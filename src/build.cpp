#include "covidx/detail/build.h"

#include <vector>

#include "covidx/detail/scope.h"

namespace covidx::detail {

uint16_t effective_exclude_flags(const BuildOptions& options) {
  uint16_t mask = options.exclude_flags;
  mask = static_cast<uint16_t>(mask | 0x0800);
  if (options.include_supplementary)
    mask = static_cast<uint16_t>(mask & ~0x0800);
  return mask;
}

void build_single_unit(const BamHeaderInfo& header, const BuildOptions& options,
                       uint32_t tid, RunIndex& out) {
  if (tid >= header.contigs.size())
    throw UsageError("contig ID is out of range");
  const ContigRecord& contig = header.contigs[tid];
  build_single_unit_scoped(header, options, tid, {{0, contig.length}}, out);
}

} // namespace covidx::detail
