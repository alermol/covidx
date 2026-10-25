#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include <htslib/sam.h>

#include "covidx/detail/cigar.h"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, std::size_t size) {
  if (size < 2 || size > 4096)
    return 0;
  const std::size_t cigar_bytes = (size - 2) / 4;
  if (cigar_bytes == 0)
    return 0;
  std::vector<uint32_t> cigar(cigar_bytes);
  std::memcpy(cigar.data(), data + 2, cigar.size() * 4);
  const int32_t pos = static_cast<int32_t>(data[0] | (data[1] << 8));
  const uint64_t contig_length = 1u << 20;
  if (pos < 0)
    return 0;

  bam1_t* record = bam_init1();
  if (record == nullptr)
    return 0;
  const std::string seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
  const std::string qual(seq.size(), 'I');
  const int rc =
      bam_set1(record, 4, "fuzz", 0, 0, pos, 60, cigar.size(), cigar.data(), -1,
               -1, 0, seq.size(), seq.c_str(), qual.c_str(), 0);
  if (rc >= 0) {
    try {
      (void)covidx::detail::extract_blocks(
          *record,
          covidx::detail::ExtractPolicy{true,
                                        covidx::AmbiguousBasePolicy::Cover},
          contig_length);
    } catch (...) {
    }
  }
  bam_destroy1(record);
  return 0;
}
