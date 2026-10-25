#pragma once
//
// Public covidx API: build, save, load, and O(1) coverage queries.
//
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "covidx/error.h"
#include "covidx/format.h"

namespace covidx {

class Manifest;

enum class LogLevel { Info, Warning };
using LogSink = std::function<void(LogLevel, std::string_view)>;

struct BuildOptions {
  unsigned threads = 0; // 0 = auto
  uint8_t min_mapq = 0;
  uint16_t exclude_flags = 0x704;
  bool include_supplementary = false;
  bool count_deletions = true;
  AmbiguousBasePolicy ambiguous_read_bases = AmbiguousBasePolicy::Cover;
  Layout layout = Layout::PerSample;
  std::vector<Group> groups;
  std::string unknown_sample;
  std::vector<std::string> contigs;
  std::vector<Region> regions;
  bool checksum_source = false;
  std::string work_dir;
  bool keep_work = false;
  bool verbose = false;
  bool require_fresh_index = false;
  LogSink log;
};

struct OpenOptions {
  bool validate_crc = true;
  bool mmap = true;
  Validation validation = Validation::Structural;
  LogSink log;
};

struct BuildStats {
  uint64_t records_examined = 0;
  uint64_t records_unmapped = 0;
  uint64_t records_flag_filtered = 0;
  uint64_t records_mapq_filtered = 0;
  uint64_t records_no_blocks = 0;
  uint64_t records_accepted = 0;
  uint64_t runs = 0;
  uint64_t units = 0;
  uint64_t contigs = 0;
};

// Immutable coverage index. Build it from a BAM or load a published
// directory; queries are thread-safe once build/load returns.
class CovIndex {
public:
  CovIndex();
  ~CovIndex();
  CovIndex(CovIndex&&) noexcept;
  CovIndex& operator=(CovIndex&&) noexcept;
  CovIndex(const CovIndex&) = delete;
  CovIndex& operator=(const CovIndex&) = delete;

  void build(const std::string& bam_path, const BuildOptions& options = {});
  void save(const std::string& out_dir, const std::string& name,
            bool overwrite = false) const;
  void load(const std::string& dir, const OpenOptions& options = {});

  [[nodiscard]] uint64_t covered_count(std::string_view unit,
                                       std::string_view contig, uint64_t x,
                                       uint64_t y) const;
  [[nodiscard]] double covered_fraction(std::string_view unit,
                                        std::string_view contig, uint64_t x,
                                        uint64_t y) const;
  [[nodiscard]] uint64_t total_covered(std::string_view unit,
                                       std::string_view contig) const;

  [[nodiscard]] const std::vector<UnitInfo>& units() const;
  [[nodiscard]] const std::vector<ContigInfo>& contigs() const;
  [[nodiscard]] const UnitInfo* unit(std::string_view name) const;
  [[nodiscard]] const ContigInfo* contig(std::string_view name) const;
  [[nodiscard]] const std::vector<Region>& scope(std::string_view contig) const;
  [[nodiscard]] const Manifest& manifest() const;
  [[nodiscard]] const BuildStats& build_stats() const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace covidx
