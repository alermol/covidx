#pragma once
//
// Build-scope planning and block clipping.
//
#include <cstdint>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "covidx/cov_index.h"
#include "covidx/detail/bam_reader.h"
#include "covidx/detail/cigar.h"
#include "covidx/detail/run_index.h"

namespace covidx::detail {

struct Scope {
  // Sorted by BAM header contig order, then start. Intervals of one contig are
  // merged and non-adjacent.
  std::vector<Region> intervals;
  bool partial = false;
};

// Builds the scope from --contigs/--regions. Throws UsageError for invalid
// --contigs options and mutually exclusive filters; malformed --regions data
// (unknown contig, empty/reversed/out-of-range interval) is InputError.
[[nodiscard]] Scope make_scope(const BamHeaderInfo& header,
                               const BuildOptions& options);

// Intersects sorted blocks with sorted, disjoint scope intervals. The result
// is sorted, disjoint, and non-adjacent.
[[nodiscard]] std::vector<Block>
clip_blocks(const std::vector<Block>& blocks,
            const std::vector<std::pair<uint64_t, uint64_t>>& intervals);

// True when [x, y) is fully contained in one scope interval.
[[nodiscard]] bool window_in_scope(const std::vector<Region>& intervals,
                                   uint64_t x, uint64_t y);

// Returns the merged intervals of one contig in a scope.
[[nodiscard]] std::vector<std::pair<uint64_t, uint64_t>>
scope_intervals(const Scope& scope, std::string_view contig);

// Reads one contig through the merged scope intervals and merges accepted
// blocks into a full-length run index. Positions outside the scope are never
// indexed.
void build_single_unit_scoped(
    const BamHeaderInfo& header, const BuildOptions& options, uint32_t tid,
    const std::vector<std::pair<uint64_t, uint64_t>>& intervals, RunIndex& out);

} // namespace covidx::detail
