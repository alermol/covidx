#pragma once
//
// Unit file writer/reader. One self-contained file per unit:
//   64-byte header, N * 48-byte metadata records, aligned coverage segments,
//   and optional feature-count segments.
//
#include <array>
#include <cstdint>
#include <cstdio>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "covidx/detail/bucket_table.h"
#include "covidx/detail/run.h"
#include "covidx/format.h"

namespace covidx::detail {

struct UnitHeader {
  uint32_t unit_id = 0;
  uint32_t contig_count = 0;
  uint32_t meta_record_size = kMetaRecordSize;
  uint32_t flags = 0;
  uint64_t metadata_offset = 64;
  uint64_t data_offset = 64;
  uint64_t file_size = 0;
};

struct UnitContigMeta {
  uint64_t segment_offset = 0;
  uint64_t runs = 0;
  uint64_t total_covered = 0;
  uint64_t max_bucket_population = 0;
  uint32_t w = 20;
  uint8_t flags = 0;
  // 0 when the contig has no feature segment.
  uint64_t feature_offset = 0;

  [[nodiscard]] std::array<uint8_t, kMetaRecordSize> encode() const;
  [[nodiscard]] static UnitContigMeta decode(std::span<const uint8_t> bytes);
};

[[nodiscard]] uint32_t read_crc_sidecar(const std::string& path);
void write_crc_sidecar(const std::string& path, uint32_t crc);

class UnitFileWriter {
public:
  UnitFileWriter(std::string path, uint32_t unit_id, uint32_t contig_count,
                 bool feature_counts = false);
  ~UnitFileWriter();

  UnitFileWriter(const UnitFileWriter&) = delete;
  UnitFileWriter& operator=(const UnitFileWriter&) = delete;

  void append_contig(const UnitContigMeta& meta,
                     std::span<const uint32_t> start,
                     std::span<const uint32_t> end,
                     std::span<const uint64_t> cum,
                     std::span<const uint32_t> first,
                     std::span<const uint32_t> feature_start = {},
                     std::span<const uint32_t> feature_end = {});
  void write_meta();
  [[nodiscard]] uint32_t close_and_crc();

private:
  void write_bytes(const void* data, std::size_t size);
  void align8();
  void write_feature_segment(std::span<const uint32_t> starts,
                             std::span<const uint32_t> ends);

  std::string path_;
  std::FILE* file_ = nullptr;
  uint32_t unit_id_ = 0;
  uint32_t contig_count_ = 0;
  uint64_t offset_ = 0;
  std::vector<UnitContigMeta> metadata_;
  bool feature_counts_ = false;
  bool finished_ = false;
};

class UnitFileReader {
public:
  // contig_lengths supplies the lengths used by the query guard and bucket
  // validation. If mmap is true and the platform allows it, arrays are read
  // zero-copy; otherwise the file is read into memory.
  UnitFileReader(std::string path, std::vector<uint64_t> contig_lengths,
                 bool use_mmap, Validation level = Validation::Structural);
  ~UnitFileReader();

  UnitFileReader(const UnitFileReader&) = delete;
  UnitFileReader& operator=(const UnitFileReader&) = delete;

  [[nodiscard]] const UnitHeader& header() const {
    return header_;
  }
  [[nodiscard]] const UnitContigMeta& meta(uint32_t contig_id) const;
  [[nodiscard]] RunIndexView view(uint32_t contig_id) const;
  [[nodiscard]] bool has_feature_counts() const {
    return (header_.flags & kUnitFlagFeatureCounts) != 0;
  }
  [[nodiscard]] uint64_t feature_count(uint32_t contig_id, uint64_t x,
                                       uint64_t y) const;
  void feature_counts(uint32_t contig_id,
                      std::span<const std::pair<uint64_t, uint64_t>> windows,
                      std::span<uint64_t> out, bool add) const;
  [[nodiscard]] uint64_t feature_total(uint32_t contig_id) const;
  [[nodiscard]] bool used_mmap() const {
    return used_mmap_;
  }
  [[nodiscard]] uint32_t file_crc32() const;

  void validate(Validation level) const;

private:
  struct SwappedArrays {
    std::vector<uint32_t> start;
    std::vector<uint32_t> end;
    std::vector<uint64_t> cum;
    std::vector<uint32_t> first;
  };

  struct FeatureData {
    uint64_t offset = 0;
    uint64_t end = 0;
    uint64_t count = 0;
    uint64_t block_size = 0;
    uint64_t block_count = 0;
    uint64_t starts_offset = 0;
    uint64_t starts_bytes = 0;
    uint64_t ends_offset = 0;
    uint64_t ends_bytes = 0;
  };

  [[nodiscard]] const uint8_t* bytes_at(uint64_t offset, uint64_t size) const;
  void load_arrays(uint32_t contig_id, RunIndexView& view) const;
  void validate_contig(uint32_t contig_id, Validation level) const;
  void parse_features(uint32_t contig_id);
  void validate_features(uint32_t contig_id, Validation level) const;
  [[nodiscard]] uint64_t feature_rank(uint32_t contig_id, bool starts,
                                      uint64_t value, bool upper) const;
  [[nodiscard]] uint64_t
  decode_feature_block(const FeatureData& data, bool starts, uint64_t block,
                       std::array<uint32_t, kFeatureBlockSize>& values) const;
  [[nodiscard]] static uint32_t decode_varint(const uint8_t*& cursor,
                                              const uint8_t* end);
  [[nodiscard]] uint64_t feature_block_offset(const FeatureData& data,
                                              bool starts,
                                              uint64_t block) const;
  [[nodiscard]] uint32_t feature_block_first(const FeatureData& data,
                                             bool starts, uint64_t block) const;
  [[nodiscard]] uint64_t feature_section_offset(const FeatureData& data,
                                                bool starts) const;
  [[nodiscard]] uint64_t feature_section_bytes(const FeatureData& data,
                                               bool starts) const;

  std::string path_;
  std::vector<uint64_t> contig_lengths_;
  uint64_t file_size_ = 0;
  UnitHeader header_;
  std::vector<UnitContigMeta> metadata_;
  std::vector<FeatureData> features_;
  void* mapping_ = nullptr;
  std::vector<uint8_t> heap_;
  const uint8_t* data_ = nullptr;
  bool used_mmap_ = false;
  mutable std::vector<SwappedArrays> swapped_;
};

} // namespace covidx::detail
