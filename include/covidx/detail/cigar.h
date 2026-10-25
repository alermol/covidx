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

struct ReferenceSpan {
  uint64_t start = 0;
  uint64_t end = 0;
};

// Returns sorted, non-overlapping, non-adjacent covered blocks within
// [0, contig_length). When span is non-null, stores the full reference span
// consumed by M/D/N/=/X. Throws InputError for malformed CIGAR, out-of-range
// coordinates, or SEQ=* when the policy needs the sequence.
[[nodiscard]] std::vector<Block> extract_blocks(const bam1_t& read,
                                                const ExtractPolicy& policy,
                                                uint64_t contig_length,
                                                ReferenceSpan* span = nullptr);

struct ReadFilter {
  uint8_t min_mapq = 0;
  uint16_t effective_exclude = 0x0F04;
  ExtractPolicy extract;
  uint64_t contig_length = 0;

  // Always skips unmapped/tid<0/pos<0 records. Returns true if accepted and
  // fills blocks; throws InputError for malformed CIGAR/SEQ cases.
  [[nodiscard]] bool accept(const bam1_t& read,
                            std::vector<Block>& blocks) const;

  // Like accept(), but returns true even when the policy leaves no covered
  // blocks. The full reference span is always filled for accepted records.
  [[nodiscard]] bool accept_alignment(const bam1_t& read,
                                      std::vector<Block>& blocks,
                                      ReferenceSpan& span) const;
};

} // namespace covidx::detail
