#pragma once
//
// Internal 16-byte run spools and the owned build workspace.
//
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "covidx/cov_index.h"
#include "covidx/detail/bam_reader.h"
#include "covidx/detail/merge.h"
#include "covidx/detail/run.h"

namespace covidx::detail {

class RunSpoolWriter {
public:
  explicit RunSpoolWriter(std::string path);
  ~RunSpoolWriter();

  RunSpoolWriter(const RunSpoolWriter&) = delete;
  RunSpoolWriter& operator=(const RunSpoolWriter&) = delete;

  void write(const Run& run);
  void close();
  [[nodiscard]] uint64_t records() const {
    return records_;
  }
  [[nodiscard]] const std::string& path() const {
    return path_;
  }

private:
  void flush();

  std::string path_;
  std::vector<uint8_t> buffer_;
  std::FILE* file_ = nullptr;
  uint64_t records_ = 0;
};

// Reads a spool sequentially. A missing file is treated as an empty spool.
// Records are validated as sorted, non-overlapping, non-adjacent runs.
class RunSpoolReader {
public:
  explicit RunSpoolReader(std::string path);
  ~RunSpoolReader();

  RunSpoolReader(const RunSpoolReader&) = delete;
  RunSpoolReader& operator=(const RunSpoolReader&) = delete;

  [[nodiscard]] bool next(Run& run);
  [[nodiscard]] std::vector<Run> read_all();
  [[nodiscard]] const std::string& path() const {
    return path_;
  }

private:
  std::string path_;
  std::FILE* file_ = nullptr;
  bool have_previous_ = false;
  Run previous_;
};

// Streams the union of several sorted spools into one writer.
void merge_run_spools(std::vector<std::unique_ptr<RunSpoolReader>>& readers,
                      RunSpoolWriter& writer);

class BuildWorkspace {
public:
  explicit BuildWorkspace(std::string parent, bool keep);
  ~BuildWorkspace();

  BuildWorkspace(const BuildWorkspace&) = delete;
  BuildWorkspace& operator=(const BuildWorkspace&) = delete;

  [[nodiscard]] const std::filesystem::path& root() const {
    return root_;
  }
  [[nodiscard]] std::filesystem::path data_dir() const {
    return root_ / "data";
  }
  [[nodiscard]] std::filesystem::path spool_dir() const {
    return root_ / "spools";
  }
  [[nodiscard]] std::filesystem::path sample_spool(uint32_t sample_id,
                                                   uint32_t contig_id) const;
  [[nodiscard]] std::filesystem::path unit_spool(uint32_t unit_id,
                                                 uint32_t contig_id) const;
  [[nodiscard]] bool keep() const {
    return keep_;
  }
  [[nodiscard]] const std::string& path_string() const {
    return root_string_;
  }

  void remove_spools();

private:
  std::filesystem::path root_;
  std::string root_string_;
  bool keep_ = false;
};

// Returns the spool containing a unit's runs: the sample spool for sample
// units, or the merged unit spool for groups/pooled.
[[nodiscard]] std::filesystem::path
spool_path_for_unit(const BuildWorkspace& workspace, const UnitPlan& unit,
                    uint32_t contig_id);

// Builds coverage spools for every contig/unit in the selected scope. The
// sample spools are written by the coverage pass; group and pooled spools are
// produced by streaming merges. Deterministic for any thread count.
[[nodiscard]] BuildStats build_run_spools(
    const BamHeaderInfo& header, const BuildOptions& options,
    const std::vector<UnitPlan>& units, const std::vector<uint32_t>& bam_tids,
    const std::vector<std::vector<std::pair<uint64_t, uint64_t>>>& intervals,
    BuildWorkspace& workspace, std::size_t max_open_spools = 0);

} // namespace covidx::detail
