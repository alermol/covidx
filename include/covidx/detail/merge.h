#pragma once
//
// Unit planning and run-list union.
//
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "covidx/cov_index.h"
#include "covidx/detail/bam_reader.h"
#include "covidx/detail/run.h"

namespace covidx::detail {

struct UnitPlan {
  uint32_t id = 0;
  std::string name;
  UnitInfo::Kind kind = UnitInfo::Kind::Sample;
  std::vector<std::string> members;
  std::vector<uint32_t> sample_ids;
  bool direct_from_reads = false;
};

// Validates layout/groups and assigns deterministic unit IDs: samples sorted
// by name, then groups sorted by name, then pooled.
[[nodiscard]] std::vector<UnitPlan> make_units(const BamHeaderInfo& header,
                                               const BuildOptions& options);

// Exact union of sorted run lists. Output is canonical.
[[nodiscard]] std::vector<Run>
merge_run_lists(const std::vector<const std::vector<Run>*>& lists);

// Reads one contig once and produces runs for every unit, indexed by unit ID.
// Positions outside the scope intervals are never included.
void build_contig_all_units(
    const BamHeaderInfo& header, const BuildOptions& options,
    const std::vector<UnitPlan>& units, uint32_t tid,
    const std::vector<std::pair<uint64_t, uint64_t>>& intervals,
    std::vector<std::vector<Run>>& out);

} // namespace covidx::detail
