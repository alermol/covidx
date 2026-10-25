#pragma once
//
// Unit file writer/reader. One self-contained file per unit:
//   64-byte header, N * 48-byte metadata records, then aligned segments.
//
#include <array>
#include <cstdint>
#include <cstdio>
#include <span>
#include <string>
#include <string_view>
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

  [[nodiscard]] std::array<uint8_t, kMetaRecordSize> encode() const;
  [[nodiscard]] static UnitContigMeta decode(std::span<const uint8_t> bytes);
};

[[nodiscard]] uint32_t read_crc_sidecar(const std::string& path);
void write_crc_sidecar(const std::string& path, uint32_t crc);

class UnitFileWriter {
public:
  UnitFileWriter(std::string path, uint32_t unit_id, uint32_t contig_count);
  ~UnitFileWriter();

  UnitFileWriter(const UnitFileWriter&) = delete;
  UnitFileWriter& operator=(const UnitFileWriter&) = delete;

  void append_contig(const UnitContigMeta& meta,
                     std::span<const uint32_t> start,
                     std::span<const uint32_t> end,
                     std::span<const uint64_t> cum,
                     std::span<const uint32_t> first);
  void write_meta();
  [[nodiscard]] uint32_t close_and_crc();
  [[nodiscard]] uint64_t size() const {
    return offset_;
  }

private:
  void write_bytes(const void* data, std::size_t size);
  void align8();

  std::string path_;
  std::FILE* file_ = nullptr;
  uint32_t unit_id_ = 0;
  uint32_t contig_count_ = 0;
  uint64_t offset_ = 0;
  std::vector<UnitContigMeta> metadata_;
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

  [[nodiscard]] const uint8_t* bytes_at(uint64_t offset, uint64_t size) const;
  void load_arrays(uint32_t contig_id, RunIndexView& view) const;
  void validate_contig(uint32_t contig_id, Validation level) const;

  std::string path_;
  std::vector<uint64_t> contig_lengths_;
  uint64_t file_size_ = 0;
  UnitHeader header_;
  std::vector<UnitContigMeta> metadata_;
  void* mapping_ = nullptr;
  std::vector<uint8_t> heap_;
  const uint8_t* data_ = nullptr;
  bool used_mmap_ = false;
  mutable std::vector<SwappedArrays> swapped_;
};

} // namespace covidx::detail
