#pragma once
//
// Internal build interfaces.
//
#include <cstdint>
#include <string>
#include <vector>

#include "covidx/cov_index.h"
#include "covidx/detail/bam_reader.h"
#include "covidx/detail/cigar.h"
#include "covidx/detail/run_index.h"

namespace covidx::detail {

// Effective flag mask after applying the implicit supplementary policy.
[[nodiscard]] uint16_t effective_exclude_flags(const BuildOptions& options);

// Reads one contig in full and merges all accepted reads into one run index.
// Used by single-sample tests and the pooled-only path.
void build_single_unit(const BamHeaderInfo& header, const BuildOptions& options,
                       uint32_t tid, RunIndex& out);

} // namespace covidx::detail
