#pragma once
//
// CIGAR processing and read filtering.
//
#include <cstdint>
#include <vector>

#include <htslib/sam.h>

#include "covidx/format.h"

namespace covidx::detail {

struct Block {
  uint64_t start = 0;
  uint64_t end = 0;
};

struct ExtractPolicy {
  bool count_deletions = true;
  AmbiguousBasePolicy ambiguous = AmbiguousBasePolicy::Cover;
};

// Returns sorted, non-overlapping, non-adjacent covered blocks within
// [0, contig_length). Throws InputError for malformed CIGAR, out-of-range
// coordinates, or SEQ=* when the policy needs the sequence.
[[nodiscard]] std::vector<Block> extract_blocks(const bam1_t& read,
                                                const ExtractPolicy& policy,
                                                uint64_t contig_length);

struct ReadFilter {
  uint8_t min_mapq = 0;
  uint16_t effective_exclude = 0x0F04;
  ExtractPolicy extract;
  uint64_t contig_length = 0;

  // Always skips unmapped/tid<0/pos<0 records. Returns true if accepted and
  // fills blocks; throws InputError for malformed CIGAR/SEQ cases.
  [[nodiscard]] bool accept(const bam1_t& read,
                            std::vector<Block>& blocks) const;
};

} // namespace covidx::detail
