#include "covidx/detail/unit_file.h"

#include <algorithm>
#include <bit>
#include <cerrno>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "covidx/detail/crc32.h"
#include "covidx/error.h"

namespace covidx::detail {

namespace {

constexpr std::size_t kHeaderBytes = 64;
constexpr std::size_t kFeatureHeaderBytes = 64;
constexpr uint32_t kFeatureVersion = 1;
constexpr char kFeatureMagic[8] = {'C', 'O', 'V', 'I', 'D', 'X', 'F', '1'};

void put_varint(std::vector<uint8_t>& out, uint32_t value) {
  while (value >= 0x80) {
    out.push_back(static_cast<uint8_t>((value & 0x7F) | 0x80));
    value >>= 7;
  }
  out.push_back(static_cast<uint8_t>(value));
}

void put_u32(uint8_t* out, uint32_t value) {
  for (int i = 0; i < 4; ++i)
    out[i] = static_cast<uint8_t>(value >> (8 * i));
}

uint32_t get_u32(const uint8_t* in) {
  uint32_t value = 0;
  for (int i = 0; i < 4; ++i)
    value |= static_cast<uint32_t>(in[i]) << (8 * i);
  return value;
}

void put_u64(uint8_t* out, uint64_t value) {
  for (int i = 0; i < 8; ++i)
    out[i] = static_cast<uint8_t>(value >> (8 * i));
}

uint64_t get_u64(const uint8_t* in) {
  uint64_t value = 0;
  for (int i = 0; i < 8; ++i)
    value |= static_cast<uint64_t>(in[i]) << (8 * i);
  return value;
}

[[nodiscard]] uint64_t align_up8(uint64_t value) {
  return (value + 7) & ~uint64_t{7};
}

[[nodiscard]] uint64_t checked_add(uint64_t a, uint64_t b, const char* what) {
  if (a > std::numeric_limits<uint64_t>::max() - b) {
    throw FormatError(std::string("overflow while computing ") + what);
  }
  return a + b;
}

[[nodiscard]] uint64_t checked_mul(uint64_t a, uint64_t b, const char* what) {
  if (a != 0 && b > std::numeric_limits<uint64_t>::max() / a) {
    throw FormatError(std::string("overflow while computing ") + what);
  }
  return a * b;
}

} // namespace

std::array<uint8_t, kMetaRecordSize> UnitContigMeta::encode() const {
  std::array<uint8_t, kMetaRecordSize> out{};
  put_u64(out.data(), segment_offset);
  put_u64(out.data() + 8, runs);
  put_u64(out.data() + 16, total_covered);
  put_u64(out.data() + 24, max_bucket_population);
  put_u32(out.data() + 32, w);
  out[36] = flags;
  put_u64(out.data() + 40, feature_offset);
  return out;
}

UnitContigMeta UnitContigMeta::decode(std::span<const uint8_t> bytes) {
  if (bytes.size() != kMetaRecordSize) {
    throw FormatError("unit metadata record must be 48 bytes");
  }
  UnitContigMeta meta;
  meta.segment_offset = get_u64(bytes.data());
  meta.runs = get_u64(bytes.data() + 8);
  meta.total_covered = get_u64(bytes.data() + 16);
  meta.max_bucket_population = get_u64(bytes.data() + 24);
  meta.w = get_u32(bytes.data() + 32);
  meta.flags = bytes[36];
  meta.feature_offset = get_u64(bytes.data() + 40);
  return meta;
}

uint32_t read_crc_sidecar(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in)
    throw IoError("cannot open CRC sidecar '" + path + "'");
  std::string text((std::istreambuf_iterator<char>(in)),
                   std::istreambuf_iterator<char>());
  if (in.bad())
    throw IoError("read error on CRC sidecar '" + path + "'");
  if (text.size() != 9 || text[8] != '\n') {
    throw FormatError("CRC sidecar '" + path +
                      "' must be 8 hex digits and one LF");
  }
  for (std::size_t i = 0; i < 8; ++i) {
    if (!((text[i] >= '0' && text[i] <= '9') ||
          (text[i] >= 'a' && text[i] <= 'f'))) {
      throw FormatError("CRC sidecar '" + path + "' must use lowercase hex");
    }
  }
  uint32_t value = 0;
  for (std::size_t i = 0; i < 8; ++i) {
    value <<= 4;
    value |= static_cast<uint32_t>(text[i] <= '9' ? text[i] - '0'
                                                  : text[i] - 'a' + 10);
  }
  return value;
}

void write_crc_sidecar(const std::string& path, uint32_t crc) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out)
    throw OutputError("cannot create CRC sidecar '" + path + "'");
  const std::string text = crc32_hex(crc) + "\n";
  out.write(text.data(), static_cast<std::streamsize>(text.size()));
  if (!out)
    throw OutputError("write error on CRC sidecar '" + path + "'");
  out.close();
  if (!out)
    throw OutputError("close error on CRC sidecar '" + path + "'");
}

UnitFileWriter::UnitFileWriter(std::string path, uint32_t unit_id,
                               uint32_t contig_count, bool feature_counts)
    : path_(std::move(path)), unit_id_(unit_id), contig_count_(contig_count),
      feature_counts_(feature_counts) {
  file_ = std::fopen(path_.c_str(), "w+b");
  if (file_ == nullptr) {
    throw OutputError("cannot create unit file '" + path_ +
                      "': " + std::strerror(errno));
  }
  try {
    const uint64_t meta_bytes =
        checked_mul(contig_count_, kMetaRecordSize, "metadata size");
    const uint64_t data_offset =
        align_up8(checked_add(kHeaderBytes, meta_bytes, "data offset"));
    uint8_t header[kHeaderBytes]{};
    std::memcpy(header, kUnitMagic, 8);
    put_u32(header + 8, unit_id_);
    put_u32(header + 12, contig_count_);
    put_u32(header + 16, kMetaRecordSize);
    put_u32(header + 20, feature_counts_ ? kUnitFlagFeatureCounts : 0);
    put_u64(header + 24, kHeaderBytes);
    put_u64(header + 32, data_offset);
    write_bytes(header, sizeof(header));
    std::vector<uint8_t> zeros(
        static_cast<std::size_t>(data_offset - kHeaderBytes), 0);
    if (!zeros.empty())
      write_bytes(zeros.data(), zeros.size());
    offset_ = data_offset;
  } catch (...) {
    std::fclose(file_);
    file_ = nullptr;
    throw;
  }
}

UnitFileWriter::~UnitFileWriter() {
  if (file_ != nullptr) {
    std::fclose(file_);
    file_ = nullptr;
  }
}

void UnitFileWriter::write_bytes(const void* data, std::size_t size) {
  if (size == 0)
    return;
  if (std::fwrite(data, 1, size, file_) != size) {
    throw OutputError("write error on unit file '" + path_ + "'");
  }
  offset_ += size;
}

void UnitFileWriter::align8() {
  static const uint8_t zeros[8] = {};
  const std::size_t pad = static_cast<std::size_t>((8 - (offset_ % 8)) % 8);
  if (pad != 0)
    write_bytes(zeros, pad);
}

void UnitFileWriter::append_contig(const UnitContigMeta& input,
                                   std::span<const uint32_t> start,
                                   std::span<const uint32_t> end,
                                   std::span<const uint64_t> cum,
                                   std::span<const uint32_t> first,
                                   std::span<const uint32_t> feature_start,
                                   std::span<const uint32_t> feature_end) {
  if (metadata_.size() >= contig_count_) {
    throw OutputError("too many contig segments for unit file '" + path_ + "'");
  }
  UnitContigMeta meta = input;
  if (start.size() != end.size() || start.size() != meta.runs) {
    throw UsageError("unit file run array sizes disagree");
  }
  if (feature_start.size() != feature_end.size()) {
    throw UsageError("unit file feature array sizes disagree");
  }
  if (!feature_start.empty() && !feature_counts_) {
    throw UsageError("unit file was not opened for feature counts");
  }
  for (std::size_t i = 0; i < feature_start.size(); ++i) {
    if (i > 0 && feature_start[i - 1] > feature_start[i]) {
      throw UsageError("unit file feature starts are not sorted");
    }
    if (i > 0 && feature_end[i - 1] > feature_end[i]) {
      throw UsageError("unit file feature ends are not sorted");
    }
    if (feature_start[i] >= feature_end[i]) {
      throw UsageError("unit file feature span is empty or reversed");
    }
  }
  if (meta.runs == 0) {
    if (!start.empty() || !end.empty() || !first.empty()) {
      throw UsageError("zero-run contig has arrays");
    }
    if (!cum.empty() && (cum.size() != 1 || cum[0] != 0)) {
      throw UsageError("zero-run contig has a non-trivial cum array");
    }
    meta.segment_offset = 0;
    meta.total_covered = 0;
    meta.max_bucket_population = 0;
    meta.w = 20;
    meta.feature_offset = 0;
    if (!feature_start.empty()) {
      align8();
      meta.feature_offset = offset_;
      write_feature_segment(feature_start, feature_end);
      align8();
    }
    metadata_.push_back(meta);
    return;
  }
  if (cum.size() != start.size() + 1)
    throw UsageError("cum array must have R+1 entries");
  if (first.empty())
    throw UsageError("non-empty run set has an empty bucket table");
  uint64_t total = 0;
  for (std::size_t i = 0; i < start.size(); ++i) {
    if (start[i] >= end[i])
      throw UsageError("unit file run is empty or reversed");
    if (i > 0 && end[i - 1] >= start[i])
      throw UsageError("unit file runs overlap or touch");
    total += end[i] - start[i];
  }
  if (cum.front() != 0 || cum.back() != total || total != meta.total_covered) {
    throw UsageError("unit file cum array is inconsistent");
  }
  align8();
  meta.segment_offset = offset_;
  std::vector<uint8_t> buffer;
  buffer.resize(start.size() * 4 + end.size() * 4 + cum.size() * 8 +
                first.size() * 4);
  std::size_t at = 0;
  for (uint32_t value : start) {
    put_u32(buffer.data() + at, value);
    at += 4;
  }
  for (uint32_t value : end) {
    put_u32(buffer.data() + at, value);
    at += 4;
  }
  for (uint64_t value : cum) {
    put_u64(buffer.data() + at, value);
    at += 8;
  }
  for (uint32_t value : first) {
    put_u32(buffer.data() + at, value);
    at += 4;
  }
  write_bytes(buffer.data(), buffer.size());
  align8();
  meta.feature_offset = 0;
  if (!feature_start.empty()) {
    meta.feature_offset = offset_;
    write_feature_segment(feature_start, feature_end);
    align8();
  }
  metadata_.push_back(meta);
}

void UnitFileWriter::write_feature_segment(std::span<const uint32_t> starts,
                                           std::span<const uint32_t> ends) {
  if (starts.size() != ends.size() || starts.empty())
    throw UsageError("feature segment requires equal non-empty arrays");
  const uint64_t count = starts.size();
  const uint64_t segment_start = offset_;
  const uint64_t blocks = (count + kFeatureBlockSize - 1) / kFeatureBlockSize;

  const auto encode_section = [blocks](std::span<const uint32_t> values,
                                       uint64_t section_offset) {
    std::vector<uint64_t> offsets(static_cast<std::size_t>(blocks), 0);
    std::vector<uint32_t> lasts(static_cast<std::size_t>(blocks), 0);
    std::vector<uint8_t> encoded;
    encoded.reserve(values.size() * 2);
    const uint64_t data_start = section_offset + blocks * 8 + blocks * 4;
    for (uint64_t block = 0; block < blocks; ++block) {
      offsets[static_cast<std::size_t>(block)] = data_start + encoded.size();
      const uint64_t begin = block * kFeatureBlockSize;
      const uint64_t end =
          std::min<uint64_t>(values.size(), begin + kFeatureBlockSize);
      lasts[static_cast<std::size_t>(block)] = values[end - 1];
      uint32_t previous = 0;
      for (uint64_t i = begin; i < end; ++i) {
        const uint32_t value = values[i];
        if (i == begin) {
          put_varint(encoded, value);
        } else {
          put_varint(encoded, value - previous);
        }
        previous = value;
      }
    }
    std::vector<uint8_t> section;
    section.reserve(
        static_cast<std::size_t>(blocks * 8 + blocks * 4 + encoded.size()));
    for (uint64_t offset : offsets) {
      const std::size_t at = section.size();
      section.resize(at + 8);
      put_u64(section.data() + at, offset);
    }
    for (uint32_t last : lasts) {
      const std::size_t at = section.size();
      section.resize(at + 4);
      put_u32(section.data() + at, last);
    }
    section.insert(section.end(), encoded.begin(), encoded.end());
    return section;
  };

  const uint64_t starts_offset = kFeatureHeaderBytes;
  const std::vector<uint8_t> start_section =
      encode_section(starts, starts_offset);
  const uint64_t ends_offset = align_up8(starts_offset + start_section.size());
  const std::vector<uint8_t> end_section = encode_section(ends, ends_offset);
  const uint64_t ends_bytes = end_section.size();

  uint8_t header[kFeatureHeaderBytes]{};
  std::memcpy(header, kFeatureMagic, 8);
  put_u32(header + 8, kFeatureVersion);
  put_u32(header + 12, kFeatureBlockSize);
  put_u64(header + 16, count);
  put_u64(header + 24, starts_offset);
  put_u64(header + 32, ends_offset);
  put_u64(header + 40, start_section.size());
  put_u64(header + 48, ends_bytes);
  write_bytes(header, sizeof(header));
  write_bytes(start_section.data(), start_section.size());
  static const uint8_t zeros[8] = {};
  const std::size_t pad = static_cast<std::size_t>((8 - (offset_ % 8)) % 8);
  if (pad != 0)
    write_bytes(zeros, pad);
  if (offset_ != segment_start + ends_offset) {
    throw OutputError("feature segment alignment is inconsistent");
  }
  write_bytes(end_section.data(), end_section.size());
}

void UnitFileWriter::write_meta() {
  if (metadata_.size() != contig_count_) {
    throw OutputError("unit file '" + path_ + "' is missing contig metadata");
  }
  if (std::fseek(file_, static_cast<long>(kHeaderBytes), SEEK_SET) != 0) {
    throw OutputError("seek error on unit file '" + path_ + "'");
  }
  for (const UnitContigMeta& meta : metadata_) {
    const auto bytes = meta.encode();
    if (std::fwrite(bytes.data(), 1, bytes.size(), file_) != bytes.size()) {
      throw OutputError("write error on unit metadata '" + path_ + "'");
    }
  }
  if (std::fseek(file_, 0, SEEK_END) != 0) {
    throw OutputError("seek error on unit file '" + path_ + "'");
  }
  const long end = std::ftell(file_);
  if (end < 0)
    throw OutputError("tell error on unit file '" + path_ + "'");
  offset_ = static_cast<uint64_t>(end);
}

uint32_t UnitFileWriter::close_and_crc() {
  if (finished_)
    return crc32_file(path_);
  write_meta();
  if (std::fflush(file_) != 0) {
    throw OutputError("flush error on unit file '" + path_ + "'");
  }
  if (::fsync(::fileno(file_)) != 0) {
    throw OutputError("fsync error on unit file '" + path_ + "'");
  }
  if (std::fclose(file_) != 0) {
    file_ = nullptr;
    throw OutputError("close error on unit file '" + path_ + "'");
  }
  file_ = nullptr;
  finished_ = true;
  const uint32_t crc = crc32_file(path_);
  write_crc_sidecar(path_ + ".crc32", crc);
  return crc;
}

UnitFileReader::UnitFileReader(std::string path,
                               std::vector<uint64_t> contig_lengths,
                               bool use_mmap, Validation level)
    : path_(std::move(path)), contig_lengths_(std::move(contig_lengths)) {
  struct stat st = {};
  if (::stat(path_.c_str(), &st) != 0) {
    throw IoError("cannot stat unit file '" + path_ + "'");
  }
  file_size_ = static_cast<uint64_t>(st.st_size);
  if (file_size_ < kHeaderBytes)
    throw FormatError("unit file is shorter than its header");

  try {
    if (use_mmap) {
      const int fd = ::open(path_.c_str(), O_RDONLY);
      if (fd >= 0) {
        void* mapping = ::mmap(nullptr, static_cast<std::size_t>(file_size_),
                               PROT_READ, MAP_PRIVATE, fd, 0);
        if (mapping != MAP_FAILED) {
          mapping_ = mapping;
          data_ = static_cast<const uint8_t*>(mapping);
          used_mmap_ = true;
        }
        ::close(fd);
      }
    }
    if (data_ == nullptr) {
      std::ifstream in(path_, std::ios::binary);
      if (!in)
        throw IoError("cannot open unit file '" + path_ + "'");
      heap_.assign(std::istreambuf_iterator<char>(in),
                   std::istreambuf_iterator<char>());
      if (in.bad())
        throw IoError("read error on unit file '" + path_ + "'");
      if (heap_.size() != file_size_)
        throw FormatError("unit file size changed while reading");
      data_ = heap_.data();
    }

    if (std::memcmp(data_, kUnitMagic, 8) != 0) {
      throw FormatError("bad unit file magic in '" + path_ + "'");
    }
    header_.unit_id = get_u32(data_ + 8);
    header_.contig_count = get_u32(data_ + 12);
    header_.meta_record_size = get_u32(data_ + 16);
    header_.flags = get_u32(data_ + 20);
    header_.metadata_offset = get_u64(data_ + 24);
    header_.data_offset = get_u64(data_ + 32);
    header_.file_size = file_size_;
    if (header_.meta_record_size != kMetaRecordSize) {
      throw FormatError("unsupported unit metadata record size");
    }
    if ((header_.flags & ~kUnitFlagFeatureCounts) != 0)
      throw FormatError("unsupported unit file flags");
    if (header_.contig_count != contig_lengths_.size()) {
      throw FormatError("unit file contig count disagrees with the manifest");
    }
    const uint64_t expected_meta = kHeaderBytes;
    const uint64_t meta_bytes =
        checked_mul(header_.contig_count, kMetaRecordSize, "metadata size");
    const uint64_t expected_data =
        align_up8(checked_add(expected_meta, meta_bytes, "data offset"));
    if (header_.metadata_offset != expected_meta ||
        header_.data_offset != expected_data) {
      throw FormatError("unit file header offsets are inconsistent");
    }
    if (header_.data_offset > file_size_)
      throw FormatError("unit file data offset is out of range");

    metadata_.reserve(header_.contig_count);
    for (uint32_t i = 0; i < header_.contig_count; ++i) {
      const uint64_t offset =
          kHeaderBytes + static_cast<uint64_t>(i) * kMetaRecordSize;
      metadata_.push_back(UnitContigMeta::decode(std::span<const uint8_t>(
          bytes_at(offset, kMetaRecordSize), kMetaRecordSize)));
    }
    features_.resize(header_.contig_count);
    for (uint32_t contig = 0; contig < header_.contig_count; ++contig)
      parse_features(contig);
    swapped_.resize(header_.contig_count);
    validate(level);
    if constexpr (std::endian::native != std::endian::little) {
      // Decode the byte-swapped arrays once, while the object is still owned
      // by the constructing thread. view() then only reads them, preserving
      // the const/thread-safe query contract. Structural/Full validation has
      // already populated every non-empty contig; this also covers None.
      for (uint32_t contig = 0; contig < header_.contig_count; ++contig) {
        if (metadata_[contig].runs != 0) {
          RunIndexView scratch;
          load_arrays(contig, scratch);
        }
      }
    }
  } catch (...) {
    // A throwing constructor does not run ~UnitFileReader; release a mapping
    // made earlier in this constructor before propagating the error.
    if (mapping_ != nullptr) {
      ::munmap(mapping_, static_cast<std::size_t>(file_size_));
      mapping_ = nullptr;
      data_ = nullptr;
    }
    throw;
  }
}

UnitFileReader::~UnitFileReader() {
  if (mapping_ != nullptr)
    ::munmap(mapping_, static_cast<std::size_t>(file_size_));
}

const uint8_t* UnitFileReader::bytes_at(uint64_t offset, uint64_t size) const {
  if (offset > file_size_ || size > file_size_ - offset) {
    throw FormatError("unit file '" + path_ + "' offset is out of range");
  }
  return data_ + offset;
}

const UnitContigMeta& UnitFileReader::meta(uint32_t contig_id) const {
  if (contig_id >= metadata_.size())
    throw FormatError("unit contig ID is out of range");
  return metadata_[contig_id];
}

uint64_t UnitFileReader::feature_section_offset(const FeatureData& data,
                                                bool starts) const {
  return starts ? data.starts_offset : data.ends_offset;
}

uint64_t UnitFileReader::feature_section_bytes(const FeatureData& data,
                                               bool starts) const {
  return starts ? data.starts_bytes : data.ends_bytes;
}

uint64_t UnitFileReader::feature_block_offset(const FeatureData& data,
                                              bool starts,
                                              uint64_t block) const {
  const uint64_t section = feature_section_offset(data, starts);
  const uint8_t* base = bytes_at(data.offset + section + block * 8, 8);
  if constexpr (std::endian::native == std::endian::little) {
    uint64_t value = 0;
    std::memcpy(&value, base, sizeof(value));
    return value;
  }
  return get_u64(base);
}

uint32_t UnitFileReader::feature_block_first(const FeatureData& data,
                                             bool starts,
                                             uint64_t block) const {
  const uint64_t first_offset =
      data.offset + feature_block_offset(data, starts, block);
  const uint64_t section_end = data.offset +
                               feature_section_offset(data, starts) +
                               feature_section_bytes(data, starts);
  const uint8_t* cursor = bytes_at(first_offset, section_end - first_offset);
  return decode_varint(cursor, cursor + (section_end - first_offset));
}

uint32_t UnitFileReader::decode_varint(const uint8_t*& cursor,
                                       const uint8_t* end) {
  uint64_t value = 0;
  for (int i = 0; i < 5; ++i) {
    if (cursor >= end)
      throw FormatError("truncated feature varint");
    const uint8_t byte = *cursor++;
    value |= static_cast<uint64_t>(byte & 0x7F) << (7 * i);
    if ((byte & 0x80) == 0) {
      if (value > std::numeric_limits<uint32_t>::max()) {
        throw FormatError("feature varint exceeds uint32");
      }
      return static_cast<uint32_t>(value);
    }
  }
  throw FormatError("feature varint is too long");
}

uint64_t UnitFileReader::decode_feature_block(
    const FeatureData& data, bool starts, uint64_t block,
    std::array<uint32_t, kFeatureBlockSize>& values) const {
  if (block >= data.block_count)
    throw FormatError("feature block index is out of range");
  const uint64_t begin = block * data.block_size;
  const uint64_t count =
      std::min<uint64_t>(data.block_size, data.count - begin);
  const uint64_t first_offset =
      data.offset + feature_block_offset(data, starts, block);
  const uint64_t section_end = data.offset +
                               feature_section_offset(data, starts) +
                               feature_section_bytes(data, starts);
  const uint64_t last_offset =
      block + 1 < data.block_count
          ? data.offset + feature_block_offset(data, starts, block + 1)
          : section_end;
  if (first_offset >= last_offset || last_offset > section_end) {
    throw FormatError("feature block range is invalid");
  }
  const uint8_t* cursor = bytes_at(first_offset, last_offset - first_offset);
  const uint8_t* end = cursor + (last_offset - first_offset);
  uint32_t previous = 0;
  for (uint64_t i = 0; i < count; ++i) {
    const uint32_t delta = decode_varint(cursor, end);
    uint64_t value = delta;
    if (i > 0)
      value += previous;
    if (value > std::numeric_limits<uint32_t>::max())
      throw FormatError("feature value overflows uint32");
    if (i > 0 && value < previous)
      throw FormatError("feature block is not sorted");
    values[static_cast<std::size_t>(i)] = static_cast<uint32_t>(value);
    previous = static_cast<uint32_t>(value);
  }
  if (cursor != end)
    throw FormatError("feature block has trailing bytes");
  const uint64_t lasts_offset =
      data.offset + feature_section_offset(data, starts) + data.block_count * 8;
  const uint8_t* lasts = bytes_at(lasts_offset + block * 4, 4);
  if (get_u32(lasts) != values[static_cast<std::size_t>(count - 1)]) {
    throw FormatError("feature block last value disagrees with its index");
  }
  return count;
}

uint64_t UnitFileReader::feature_rank(uint32_t contig_id, bool starts,
                                      uint64_t value, bool upper) const {
  const FeatureData& data = features_.at(contig_id);
  if (data.count == 0)
    return 0;
  const uint64_t lasts_offset =
      data.offset + feature_section_offset(data, starts) + data.block_count * 8;
  const uint8_t* lasts_bytes = bytes_at(lasts_offset, data.block_count * 4);
  uint64_t low = 0;
  uint64_t high = data.block_count;
  if constexpr (std::endian::native == std::endian::little) {
    const uint32_t* lasts = reinterpret_cast<const uint32_t*>(lasts_bytes);
    while (low < high) {
      const uint64_t mid = low + (high - low) / 2;
      const uint32_t last = lasts[mid];
      const bool before = upper ? last <= value : last < value;
      if (before)
        low = mid + 1;
      else
        high = mid;
    }
  } else {
    while (low < high) {
      const uint64_t mid = low + (high - low) / 2;
      const uint32_t last = get_u32(lasts_bytes + mid * 4);
      const bool before = upper ? last <= value : last < value;
      if (before)
        low = mid + 1;
      else
        high = mid;
    }
  }
  if (low == data.block_count)
    return data.count;
  const uint64_t block = low;
  const uint32_t first_value = feature_block_first(data, starts, block);
  if (upper ? first_value > value : first_value >= value)
    return block * data.block_size;
  std::array<uint32_t, kFeatureBlockSize> values{};
  const uint64_t count = decode_feature_block(data, starts, block, values);
  const uint64_t begin = block * data.block_size;
  const auto first = values.begin();
  const auto last = first + static_cast<std::ptrdiff_t>(count);
  const auto it = upper ? std::upper_bound(first, last, value)
                        : std::lower_bound(first, last, value);
  return begin + static_cast<uint64_t>(it - first);
}

void UnitFileReader::feature_counts(
    uint32_t contig_id, std::span<const std::pair<uint64_t, uint64_t>> windows,
    std::span<uint64_t> out, bool add) const {
  if (contig_id >= features_.size())
    throw FormatError("unit contig ID is out of range");
  if (out.size() != windows.size())
    throw UsageError("feature batch size mismatch");
  if (!has_feature_counts())
    throw FormatError("unit file has no feature counts");
  const uint64_t length = contig_lengths_.at(contig_id);
  for (const auto& [x, y] : windows) {
    if (x >= y) {
      throw QueryError("empty or reversed feature window [" +
                       std::to_string(x) + ", " + std::to_string(y) + ")");
    }
    if (y > length) {
      throw QueryError("feature window end " + std::to_string(y) +
                       " exceeds contig length " + std::to_string(length));
    }
  }
  const FeatureData& data = features_[contig_id];
  if (data.count == 0) {
    if (!add)
      std::fill(out.begin(), out.end(), uint64_t{0});
    return;
  }

  std::vector<uint32_t> starts(static_cast<std::size_t>(data.count));
  std::vector<uint32_t> ends(static_cast<std::size_t>(data.count));
  std::array<uint32_t, kFeatureBlockSize> block_values{};
  for (uint64_t block = 0; block < data.block_count; ++block) {
    const uint64_t begin = block * data.block_size;
    const std::size_t at = static_cast<std::size_t>(begin);
    const uint64_t count =
        decode_feature_block(data, true, block, block_values);
    std::copy_n(block_values.begin(), static_cast<std::size_t>(count),
                starts.begin() + static_cast<std::ptrdiff_t>(at));
    const uint64_t end_count =
        decode_feature_block(data, false, block, block_values);
    if (end_count != count)
      throw FormatError("feature block value counts disagree");
    std::copy_n(block_values.begin(), static_cast<std::size_t>(end_count),
                ends.begin() + static_cast<std::ptrdiff_t>(at));
  }

  for (std::size_t i = 0; i < windows.size(); ++i) {
    const uint64_t x = windows[i].first;
    const uint64_t y = windows[i].second;
    const uint64_t starts_before_y = static_cast<uint64_t>(
        std::lower_bound(starts.begin(), starts.end(), y) - starts.begin());
    const uint64_t ends_at_or_before_x = static_cast<uint64_t>(
        std::upper_bound(ends.begin(), ends.end(), x) - ends.begin());
    if (starts_before_y < ends_at_or_before_x)
      throw FormatError("feature arrays are inconsistent");
    const uint64_t value = starts_before_y - ends_at_or_before_x;
    if (add)
      out[i] += value;
    else
      out[i] = value;
  }
}

uint64_t UnitFileReader::feature_total(uint32_t contig_id) const {
  if (contig_id >= features_.size())
    throw FormatError("unit contig ID is out of range");
  return features_[contig_id].count;
}

uint64_t UnitFileReader::feature_count(uint32_t contig_id, uint64_t x,
                                       uint64_t y) const {
  if (contig_id >= features_.size())
    throw FormatError("unit contig ID is out of range");
  if (x >= y) {
    throw QueryError("empty or reversed feature window [" + std::to_string(x) +
                     ", " + std::to_string(y) + ")");
  }
  const uint64_t length = contig_lengths_.at(contig_id);
  if (y > length) {
    throw QueryError("feature window end " + std::to_string(y) +
                     " exceeds contig length " + std::to_string(length));
  }
  if (!has_feature_counts())
    throw FormatError("unit file has no feature counts");
  if (features_[contig_id].count == 0)
    return 0;
  const uint64_t start_rank = feature_rank(contig_id, true, y, false);
  const uint64_t end_rank = feature_rank(contig_id, false, x, true);
  if (start_rank < end_rank) {
    throw FormatError("feature arrays are inconsistent");
  }
  return start_rank - end_rank;
}

void UnitFileReader::parse_features(uint32_t contig_id) {
  if (contig_id >= features_.size())
    throw FormatError("unit contig ID is out of range");
  FeatureData& data = features_[contig_id];
  data = FeatureData{};
  const UnitContigMeta& info = meta(contig_id);
  if (info.feature_offset == 0)
    return;
  if (!has_feature_counts())
    throw FormatError("unit has a feature segment without the feature flag");
  if ((info.feature_offset % kPackAlignment) != 0)
    throw FormatError("feature segment offset is not 8-byte aligned");
  if (info.feature_offset < header_.data_offset) {
    throw FormatError("feature segment offset is before the data area");
  }
  if (info.feature_offset > file_size_ ||
      kFeatureHeaderBytes > file_size_ - info.feature_offset) {
    throw FormatError("feature segment header is out of range");
  }

  const uint8_t* header = bytes_at(info.feature_offset, kFeatureHeaderBytes);
  if (std::memcmp(header, kFeatureMagic, 8) != 0)
    throw FormatError("bad feature segment magic");
  const uint32_t version = get_u32(header + 8);
  const uint32_t block_size = get_u32(header + 12);
  if (version != kFeatureVersion)
    throw FormatError("unsupported feature segment version");
  if (block_size != kFeatureBlockSize)
    throw FormatError("unsupported feature block size");
  data.offset = info.feature_offset;
  data.count = get_u64(header + 16);
  data.starts_offset = get_u64(header + 24);
  data.ends_offset = get_u64(header + 32);
  data.starts_bytes = get_u64(header + 40);
  data.ends_bytes = get_u64(header + 48);
  if (get_u64(header + 56) != 0)
    throw FormatError("feature segment reserved field is non-zero");
  if (data.count == 0)
    throw FormatError("feature segment contains no features");
  if (data.starts_offset != kFeatureHeaderBytes)
    throw FormatError("feature starts offset is invalid");
  data.block_size = block_size;
  data.block_count = 1 + (data.count - 1) / data.block_size;
  if (data.block_count > data.starts_bytes / 12 ||
      data.block_count > data.ends_bytes / 12) {
    throw FormatError("feature block counts exceed the section sizes");
  }
  if (data.starts_bytes >
      std::numeric_limits<uint64_t>::max() - data.starts_offset) {
    throw FormatError("feature starts section overflows");
  }
  const uint64_t starts_end = data.starts_offset + data.starts_bytes;
  if (starts_end > std::numeric_limits<uint64_t>::max() - 7) {
    throw FormatError("feature starts section overflows");
  }
  const uint64_t expected_ends = align_up8(starts_end);
  if (data.ends_offset != expected_ends)
    throw FormatError("feature ends offset is inconsistent");
  if (data.ends_bytes >
      std::numeric_limits<uint64_t>::max() - data.ends_offset) {
    throw FormatError("feature ends section overflows");
  }
  const uint64_t ends_end = data.ends_offset + data.ends_bytes;
  if (ends_end > std::numeric_limits<uint64_t>::max() - 7) {
    throw FormatError("feature ends section overflows");
  }
  data.end = align_up8(ends_end);
  if (data.end > file_size_ - data.offset)
    throw FormatError("feature segment is out of range");

  const auto check_section = [&](bool starts) {
    const uint64_t section_offset = feature_section_offset(data, starts);
    const uint64_t section_bytes = feature_section_bytes(data, starts);
    const uint64_t header_bytes = data.block_count * 12;
    if (section_bytes < header_bytes)
      throw FormatError("feature section is too short");
    const uint64_t data_begin = section_offset + header_bytes;
    const uint64_t data_end = section_offset + section_bytes;
    (void)bytes_at(data.offset + section_offset, section_bytes);
    const uint8_t* offsets =
        bytes_at(data.offset + section_offset, data.block_count * 8);
    const uint8_t* lasts =
        bytes_at(data.offset + section_offset + data.block_count * 8,
                 data.block_count * 4);
    uint64_t previous_offset = 0;
    uint32_t previous_last = 0;
    for (uint64_t block = 0; block < data.block_count; ++block) {
      const uint64_t offset = get_u64(offsets + block * 8);
      if (offset < data_begin || offset >= data_end)
        throw FormatError("feature block offset is out of range");
      if (block == 0 && offset != data_begin)
        throw FormatError("first feature block offset is invalid");
      if (block > 0 && offset <= previous_offset)
        throw FormatError("feature block offsets are not increasing");
      previous_offset = offset;
      const uint32_t last = get_u32(lasts + block * 4);
      if (block > 0 && last < previous_last)
        throw FormatError("feature block last values are not sorted");
      previous_last = last;
    }
  };
  check_section(true);
  check_section(false);
}

void UnitFileReader::load_arrays(uint32_t contig_id, RunIndexView& view) const {
  const UnitContigMeta& info = meta(contig_id);
  const uint64_t runs = info.runs;
  if (std::endian::native == std::endian::little) {
    const uint8_t* base =
        bytes_at(info.segment_offset,
                 runs * 4 + runs * 4 + (runs + 1) * 8 +
                     bucket_entries(contig_lengths_.at(contig_id), info.w) * 4);
    view.start =
        std::span<const uint32_t>(reinterpret_cast<const uint32_t*>(base),
                                  static_cast<std::size_t>(runs));
    view.end = std::span<const uint32_t>(
        reinterpret_cast<const uint32_t*>(base + runs * 4),
        static_cast<std::size_t>(runs));
    view.cum = std::span<const uint64_t>(
        reinterpret_cast<const uint64_t*>(base + runs * 8),
        static_cast<std::size_t>(runs + 1));
    const uint64_t first_count =
        bucket_entries(contig_lengths_.at(contig_id), info.w);
    view.first = std::span<const uint32_t>(
        reinterpret_cast<const uint32_t*>(base + runs * 8 + (runs + 1) * 8),
        static_cast<std::size_t>(first_count));
    return;
  }

  SwappedArrays& arrays = swapped_.at(contig_id);
  if (!arrays.start.empty()) {
    view.start = arrays.start;
    view.end = arrays.end;
    view.cum = arrays.cum;
    view.first = arrays.first;
    return;
  }
  const uint8_t* base =
      bytes_at(info.segment_offset,
               runs * 4 + runs * 4 + (runs + 1) * 8 +
                   bucket_entries(contig_lengths_.at(contig_id), info.w) * 4);
  arrays.start.resize(static_cast<std::size_t>(runs));
  arrays.end.resize(static_cast<std::size_t>(runs));
  arrays.cum.resize(static_cast<std::size_t>(runs + 1));
  const uint64_t first_count =
      bucket_entries(contig_lengths_.at(contig_id), info.w);
  arrays.first.resize(static_cast<std::size_t>(first_count));
  for (uint64_t i = 0; i < runs; ++i) {
    arrays.start[static_cast<std::size_t>(i)] = get_u32(base + i * 4);
    arrays.end[static_cast<std::size_t>(i)] = get_u32(base + runs * 4 + i * 4);
  }
  for (uint64_t i = 0; i <= runs; ++i) {
    arrays.cum[static_cast<std::size_t>(i)] = get_u64(base + runs * 8 + i * 8);
  }
  for (uint64_t i = 0; i < first_count; ++i) {
    arrays.first[static_cast<std::size_t>(i)] =
        get_u32(base + runs * 8 + (runs + 1) * 8 + i * 4);
  }
  view.start = arrays.start;
  view.end = arrays.end;
  view.cum = arrays.cum;
  view.first = arrays.first;
}

RunIndexView UnitFileReader::view(uint32_t contig_id) const {
  if (contig_id >= metadata_.size())
    throw FormatError("unit contig ID is out of range");
  const UnitContigMeta& info = meta(contig_id);
  RunIndexView view;
  view.length = contig_lengths_.at(contig_id);
  view.runs = info.runs;
  view.total_covered = info.total_covered;
  view.w = info.w;
  if (info.runs == 0)
    return view;
  load_arrays(contig_id, view);
  return view;
}

void UnitFileReader::validate_contig(uint32_t contig_id,
                                     Validation level) const {
  const UnitContigMeta& info = meta(contig_id);
  const uint64_t length = contig_lengths_.at(contig_id);
  if (length > kMaxContigLength)
    throw FormatError("contig length exceeds INT32_MAX");
  if (info.total_covered > length)
    throw FormatError("unit total_covered exceeds contig length");
  if (info.runs == 0) {
    if (info.segment_offset != 0 || info.total_covered != 0 ||
        info.max_bucket_population != 0 || info.w != 20) {
      throw FormatError("zero-run unit metadata is inconsistent");
    }
    return;
  }
  if (info.w < 6 || info.w > 20)
    throw FormatError("unit bucket width is out of range");
  if ((info.segment_offset % kPackAlignment) != 0) {
    throw FormatError("unit segment offset is not 8-byte aligned");
  }
  if (info.runs > (length + 1) / 2)
    throw FormatError("unit run count is too large");
  const uint64_t first_count = bucket_entries(length, info.w);
  const uint64_t bytes = checked_mul(info.runs, 4, "start array") +
                         checked_mul(info.runs, 4, "end array") +
                         checked_mul(info.runs + 1, 8, "cum array") +
                         checked_mul(first_count, 4, "bucket array");
  if (info.segment_offset > file_size_ ||
      bytes > file_size_ - info.segment_offset) {
    throw FormatError("unit segment is out of range");
  }
  if (level == Validation::None)
    return;

  const RunIndexView data = view(contig_id);
  if (data.start.size() != info.runs || data.end.size() != info.runs ||
      data.cum.size() != info.runs + 1 || data.first.size() != first_count) {
    throw FormatError("unit array sizes disagree with metadata");
  }
  if (data.cum.front() != 0 || data.cum.back() != info.total_covered) {
    throw FormatError("unit cum endpoints disagree with metadata");
  }
  if (data.first.front() != 0 || data.first.back() != info.runs) {
    throw FormatError("unit bucket endpoints are wrong");
  }
  if (data.start.front() >= data.end.front() ||
      data.start.back() >= data.end.back() || data.end.back() > length) {
    throw FormatError("unit run endpoints are out of range");
  }

  const auto check_index = [&](std::size_t i) {
    if (data.start[i] >= data.end[i] || data.end[i] > length) {
      throw FormatError("unit run is out of range");
    }
    if (i > 0 && data.end[i - 1] >= data.start[i]) {
      throw FormatError("unit runs overlap or touch");
    }
    const uint64_t len = data.end[i] - data.start[i];
    if (data.cum[i + 1] != data.cum[i] + len)
      throw FormatError("unit cum recurrence is broken");
  };
  const auto check_bucket = [&](std::size_t b) {
    if (data.first[b] > data.first[b + 1])
      throw FormatError("unit bucket table is not monotone");
  };

  if (level == Validation::Structural) {
    for (std::size_t i :
         {std::size_t{0}, static_cast<std::size_t>(info.runs - 1)}) {
      check_index(i);
    }
    for (std::size_t b :
         {std::size_t{0}, static_cast<std::size_t>(first_count - 2),
          static_cast<std::size_t>(first_count - 1)}) {
      if (b + 1 < data.first.size())
        check_bucket(b);
    }
    if (info.runs > 1) {
      const std::size_t mid = static_cast<std::size_t>(info.runs / 2);
      check_index(mid - 1);
      check_index(mid);
    }
    return;
  }

  uint64_t cumulative = 0;
  for (std::size_t i = 0; i < data.start.size(); ++i) {
    check_index(i);
    cumulative += data.end[i] - data.start[i];
  }
  if (cumulative != info.total_covered)
    throw FormatError("unit total coverage is wrong");
  uint32_t expected = 0;
  for (std::size_t b = 0; b < data.first.size(); ++b) {
    const uint64_t limit = static_cast<uint64_t>(b) << info.w;
    while (expected < data.start.size() && data.start[expected] < limit)
      ++expected;
    if (data.first[b] != expected)
      throw FormatError("unit bucket table value is wrong");
    if (b > 0 && data.first[b] < data.first[b - 1]) {
      throw FormatError("unit bucket table is not monotone");
    }
  }
  if (expected != data.start.size())
    throw FormatError("unit bucket table misses runs");
}

void UnitFileReader::validate_features(uint32_t contig_id,
                                       Validation level) const {
  const FeatureData& data = features_.at(contig_id);
  const UnitContigMeta& info = meta(contig_id);
  const uint64_t length = contig_lengths_.at(contig_id);
  uint64_t coverage_end = 0;
  if (info.runs > 0) {
    const uint64_t first_count = bucket_entries(length, info.w);
    const uint64_t bytes =
        info.runs * 4 + info.runs * 4 + (info.runs + 1) * 8 + first_count * 4;
    coverage_end = info.segment_offset + bytes;
  }
  if (data.offset == 0)
    return;
  if (info.feature_offset < coverage_end)
    throw FormatError("feature segment overlaps coverage segment");
  if (data.offset + data.end > file_size_)
    throw FormatError("feature segment extends beyond the unit file");
  if (level == Validation::None)
    return;

  const auto check_values = [&](bool starts, uint64_t block,
                                uint64_t previous) {
    std::array<uint32_t, kFeatureBlockSize> values{};
    const uint64_t count = decode_feature_block(data, starts, block, values);
    for (uint64_t i = 0; i < count; ++i) {
      const uint32_t value = values[static_cast<std::size_t>(i)];
      if (starts) {
        if (value >= length)
          throw FormatError("feature start is beyond the contig");
      } else if (value > length) {
        throw FormatError("feature end is beyond the contig");
      }
      if (i == 0 && block > 0 && value < previous)
        throw FormatError("feature blocks are not sorted");
    }
    return values[static_cast<std::size_t>(count - 1)];
  };

  if (level == Validation::Structural) {
    for (bool starts : {true, false}) {
      (void)check_values(starts, 0, 0);
      if (data.block_count > 1) {
        const uint64_t previous = check_values(starts, data.block_count - 2, 0);
        (void)check_values(starts, data.block_count - 1, previous);
      }
    }
    return;
  }

  uint32_t previous_start = 0;
  uint32_t previous_end = 0;
  for (uint64_t block = 0; block < data.block_count; ++block) {
    std::array<uint32_t, kFeatureBlockSize> starts{};
    std::array<uint32_t, kFeatureBlockSize> ends{};
    const uint64_t count = decode_feature_block(data, true, block, starts);
    const uint64_t end_count = decode_feature_block(data, false, block, ends);
    if (count != end_count)
      throw FormatError("feature block value counts disagree");
    for (uint64_t i = 0; i < count; ++i) {
      const uint32_t start = starts[static_cast<std::size_t>(i)];
      const uint32_t end = ends[static_cast<std::size_t>(i)];
      if (start >= length || end > length || start >= end)
        throw FormatError("feature span is out of range");
      if (block > 0 && i == 0 &&
          (start < previous_start || end < previous_end)) {
        throw FormatError("feature blocks are not sorted");
      }
      previous_start = start;
      previous_end = end;
    }
  }
}

void UnitFileReader::validate(Validation level) const {
  for (uint32_t contig = 0; contig < metadata_.size(); ++contig) {
    validate_contig(contig, level);
    validate_features(contig, level);
  }

  std::vector<std::pair<uint64_t, uint64_t>> ranges;
  ranges.reserve(metadata_.size() * 2);
  for (uint32_t contig = 0; contig < metadata_.size(); ++contig) {
    const UnitContigMeta& info = metadata_[contig];
    if (info.runs > 0) {
      const uint64_t first_count =
          bucket_entries(contig_lengths_[contig], info.w);
      const uint64_t bytes =
          info.runs * 4 + info.runs * 4 + (info.runs + 1) * 8 + first_count * 4;
      ranges.push_back({info.segment_offset, info.segment_offset + bytes});
    }
    if (features_[contig].offset != 0) {
      ranges.push_back({features_[contig].offset,
                        features_[contig].offset + features_[contig].end});
    }
  }
  std::sort(ranges.begin(), ranges.end());
  for (std::size_t i = 1; i < ranges.size(); ++i) {
    if (ranges[i - 1].second > ranges[i].first) {
      throw FormatError("unit file segments overlap");
    }
  }
}

uint32_t UnitFileReader::file_crc32() const {
  return crc32_file(path_);
}

} // namespace covidx::detail
