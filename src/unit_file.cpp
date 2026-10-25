#include "covidx/detail/unit_file.h"

#include <bit>
#include <cerrno>
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
                               uint32_t contig_count)
    : path_(std::move(path)), unit_id_(unit_id), contig_count_(contig_count) {
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
    put_u32(header + 20, 0);
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
                                   std::span<const uint32_t> first) {
  if (metadata_.size() >= contig_count_) {
    throw OutputError("too many contig segments for unit file '" + path_ + "'");
  }
  UnitContigMeta meta = input;
  if (start.size() != end.size() || start.size() != meta.runs) {
    throw UsageError("unit file run array sizes disagree");
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
  metadata_.push_back(meta);
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
    if (header_.flags != 0)
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

void UnitFileReader::validate(Validation level) const {
  for (uint32_t contig = 0; contig < metadata_.size(); ++contig) {
    validate_contig(contig, level);
  }
}

uint32_t UnitFileReader::file_crc32() const {
  return crc32_file(path_);
}

} // namespace covidx::detail
