#include "covidx/detail/cigar.h"

#include <algorithm>
#include <cstdint>
#include <string>

#include "covidx/error.h"

namespace covidx::detail {

namespace {

constexpr uint32_t kCigarM = 0;
constexpr uint32_t kCigarI = 1;
constexpr uint32_t kCigarD = 2;
constexpr uint32_t kCigarN = 3;
constexpr uint32_t kCigarS = 4;
constexpr uint32_t kCigarH = 5;
constexpr uint32_t kCigarP = 6;
constexpr uint32_t kCigarEq = 7;
constexpr uint32_t kCigarX = 8;
constexpr uint32_t kCigarB = 9;

[[nodiscard]] bool is_unambiguous_base(uint8_t code) {
  return code == 1 || code == 2 || code == 4 || code == 8;
}

[[nodiscard]] std::string read_context(const bam1_t& read) {
  const char* name = bam_get_qname(&read);
  return std::string("read '") + (name != nullptr ? name : "?") +
         "' at position " + std::to_string(read.core.pos);
}

void flush_block(std::vector<Block>& blocks, bool& active, uint64_t& start,
                 uint64_t& end) {
  if (!active)
    return;
  if (end > start)
    blocks.push_back(Block{start, end});
  active = false;
}

} // namespace

std::vector<Block> extract_blocks(const bam1_t& read,
                                  const ExtractPolicy& policy,
                                  uint64_t contig_length) {
  if (read.core.pos < 0) {
    throw InputError("mapped record has a negative position: " +
                     read_context(read));
  }
  if (static_cast<uint64_t>(read.core.pos) >= contig_length) {
    throw InputError("mapped record starts beyond contig length " +
                     std::to_string(contig_length) + ": " + read_context(read));
  }
  if (read.core.n_cigar == 0) {
    throw InputError("mapped record has no CIGAR: " + read_context(read));
  }
  if (read.core.l_qseq < 0) {
    throw InputError("mapped record has a negative sequence length: " +
                     read_context(read));
  }
  const bool needs_sequence = policy.ambiguous == AmbiguousBasePolicy::Exclude;
  if (needs_sequence && read.core.l_qseq == 0) {
    throw InputError(
        "mapped record has SEQ=* but --ambiguous-read-bases exclude needs the "
        "sequence: " +
        read_context(read));
  }

  std::vector<Block> blocks;
  blocks.reserve(read.core.n_cigar);

  uint64_t ref_pos = static_cast<uint64_t>(read.core.pos);
  uint64_t query_pos = 0;
  bool active = false;
  uint64_t cs = 0;
  uint64_t ce = 0;
  const uint32_t* cigar = bam_get_cigar(&read);
  const uint8_t* seq = read.core.l_qseq > 0 ? bam_get_seq(&read) : nullptr;

  for (uint32_t i = 0; i < read.core.n_cigar; ++i) {
    const uint32_t code = cigar[i];
    const uint32_t len = code >> 4;
    const uint32_t op = code & 0x0F;
    if (len == 0) {
      throw InputError("zero-length CIGAR operation: " + read_context(read));
    }
    if (op > kCigarB) {
      throw InputError("unsupported CIGAR operation code " +
                       std::to_string(op) + ": " + read_context(read));
    }
    switch (op) {
    case kCigarM:
    case kCigarEq:
    case kCigarX: {
      if (!needs_sequence) {
        if (!active) {
          active = true;
          cs = ref_pos;
        }
        ce = ref_pos + len;
      } else {
        if (query_pos + len > static_cast<uint64_t>(read.core.l_qseq)) {
          throw InputError("CIGAR consumes more query bases than SEQ length: " +
                           read_context(read));
        }
        for (uint32_t k = 0; k < len; ++k) {
          const uint8_t base = bam_seqi(seq, static_cast<int>(query_pos + k));
          if (is_unambiguous_base(base)) {
            if (!active) {
              active = true;
              cs = ref_pos + k;
            }
            ce = ref_pos + k + 1;
          } else {
            flush_block(blocks, active, cs, ce);
          }
        }
      }
      ref_pos += len;
      query_pos += len;
      break;
    }
    case kCigarD: {
      if (policy.count_deletions) {
        if (!active) {
          active = true;
          cs = ref_pos;
        }
        ce = ref_pos + len;
      } else {
        flush_block(blocks, active, cs, ce);
      }
      ref_pos += len;
      break;
    }
    case kCigarN: {
      flush_block(blocks, active, cs, ce);
      ref_pos += len;
      break;
    }
    case kCigarI:
    case kCigarS: {
      query_pos += len;
      break;
    }
    case kCigarH:
    case kCigarP: {
      break;
    }
    case kCigarB: {
      throw InputError("CIGAR B (back) is not supported: " +
                       read_context(read));
    }
    default: {
      throw InputError("unsupported CIGAR operation: " + read_context(read));
    }
    }
    if (ref_pos > contig_length) {
      throw InputError("alignment extends beyond contig length " +
                       std::to_string(contig_length) + ": " +
                       read_context(read));
    }
    if (query_pos > static_cast<uint64_t>(read.core.l_qseq) &&
        read.core.l_qseq > 0) {
      throw InputError("CIGAR consumes more query bases than SEQ length: " +
                       read_context(read));
    }
  }
  flush_block(blocks, active, cs, ce);

  for (const Block& block : blocks) {
    if (block.start >= block.end || block.end > contig_length) {
      throw InputError(
          "alignment block [" + std::to_string(block.start) + ", " +
          std::to_string(block.end) + ") is outside contig length " +
          std::to_string(contig_length) + ": " + read_context(read));
    }
  }
  return blocks;
}

bool ReadFilter::accept(const bam1_t& read, std::vector<Block>& blocks) const {
  blocks.clear();
  if ((read.core.flag & BAM_FUNMAP) != 0 || read.core.tid < 0 ||
      read.core.pos < 0)
    return false;
  if ((read.core.flag & effective_exclude) != 0)
    return false;
  if (read.core.qual < min_mapq)
    return false;
  blocks = extract_blocks(read, extract, contig_length);
  return !blocks.empty();
}

} // namespace covidx::detail
