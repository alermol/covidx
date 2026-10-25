#include "doctest.h"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <random>
#include <vector>

#include <htslib/sam.h>

#include "covidx/detail/bam_reader.h"
#include "covidx/detail/build.h"
#include "covidx/detail/run_index.h"
#include "test_util.h"

namespace fs = std::filesystem;
using covidx::AmbiguousBasePolicy;
using covidx::BuildOptions;
using covidx::detail::BamHeaderInfo;
using covidx::detail::BamRegionReader;
using covidx::detail::effective_exclude_flags;
using covidx::detail::read_bam_header;
using covidx::detail::RunIndex;

namespace {

std::vector<bool> naive_coverage(const BamHeaderInfo& header,
                                 const BuildOptions& options, uint32_t tid) {
  const auto& contig = header.contigs[tid];
  std::vector<bool> covered(static_cast<std::size_t>(contig.length), false);
  const uint16_t exclude = effective_exclude_flags(options);
  BamRegionReader reader(header.path, static_cast<int32_t>(tid),
                         {{0, contig.length}});
  bam1_t* record = bam_init1();
  while (reader.next(*record)) {
    if ((record->core.flag & BAM_FUNMAP) != 0 || record->core.tid < 0 ||
        record->core.pos < 0) {
      continue;
    }
    if ((record->core.flag & exclude) != 0)
      continue;
    if (record->core.qual < options.min_mapq)
      continue;
    if (options.ambiguous_read_bases == AmbiguousBasePolicy::Exclude &&
        record->core.l_qseq == 0) {
      bam_destroy1(record);
      throw covidx::InputError(
          "naive reference cannot evaluate SEQ=* under exclude");
    }
    uint64_t ref = static_cast<uint64_t>(record->core.pos);
    uint64_t query = 0;
    const uint32_t* cigar = bam_get_cigar(record);
    const uint8_t* seq =
        record->core.l_qseq > 0 ? bam_get_seq(record) : nullptr;
    for (uint32_t i = 0; i < record->core.n_cigar; ++i) {
      const uint32_t len = cigar[i] >> 4;
      const uint32_t op = cigar[i] & 0x0F;
      if (op == 0 || op == 7 || op == 8) {
        for (uint32_t k = 0; k < len; ++k) {
          bool contributes = true;
          if (options.ambiguous_read_bases == AmbiguousBasePolicy::Exclude) {
            const uint8_t base = bam_seqi(seq, static_cast<int>(query + k));
            contributes = base == 1 || base == 2 || base == 4 || base == 8;
          }
          if (contributes)
            covered[static_cast<std::size_t>(ref + k)] = true;
        }
        ref += len;
        query += len;
      } else if (op == 2) {
        if (options.count_deletions) {
          for (uint32_t k = 0; k < len; ++k)
            covered[static_cast<std::size_t>(ref + k)] = true;
        }
        ref += len;
      } else if (op == 3) {
        ref += len;
      } else if (op == 1 || op == 4) {
        query += len;
      }
    }
  }
  bam_destroy1(record);
  return covered;
}

void compare_index(const RunIndex& index, const std::vector<bool>& naive,
                   uint32_t seed) {
  REQUIRE(index.length() == naive.size());
  uint64_t total = 0;
  for (bool value : naive)
    total += value ? 1 : 0;
  CHECK(index.total_covered() == total);
  for (uint64_t x = 0; x < naive.size(); ++x) {
    uint64_t running = 0;
    const uint64_t limit = std::min<uint64_t>(naive.size(), x + 100);
    for (uint64_t y = x + 1; y <= limit; ++y) {
      running += naive[static_cast<std::size_t>(y - 1)] ? 1 : 0;
      CHECK(index.covered_count(x, y) == running);
    }
  }
  std::mt19937 rng(seed);
  for (int i = 0; i < 500; ++i) {
    const uint64_t x = rng() % naive.size();
    const uint64_t y = x + 1 + rng() % (naive.size() - x);
    uint64_t want = 0;
    for (uint64_t p = x; p < y; ++p)
      want += naive[static_cast<std::size_t>(p)] ? 1 : 0;
    CHECK(index.covered_count(x, y) == want);
  }
}

} // namespace

TEST_CASE("single-unit build matches a naive bitmap across policies") {
  const fs::path bam = fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam";
  const BamHeaderInfo header = read_bam_header(bam.string());

  std::vector<BuildOptions> policies;
  BuildOptions defaults;
  policies.push_back(defaults);

  BuildOptions no_del;
  no_del.count_deletions = false;
  policies.push_back(no_del);

  BuildOptions minq;
  minq.min_mapq = 20;
  policies.push_back(minq);

  BuildOptions all_flags;
  all_flags.exclude_flags = 0;
  policies.push_back(all_flags);

  BuildOptions include_supp;
  include_supp.include_supplementary = true;
  policies.push_back(include_supp);

  BuildOptions all_primary;
  all_primary.exclude_flags = 0x604;
  policies.push_back(all_primary);

  for (const BuildOptions& policy : policies) {
    const std::vector<bool> naive = naive_coverage(header, policy, 0);
    RunIndex index;
    covidx::detail::build_single_unit(header, policy, 0, index);
    compare_index(index, naive, 42);
    index.validate();
  }
}

TEST_CASE("ambiguous-base exclude rejects SEQ=* and works without it") {
  const fs::path bam = fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam";
  const BamHeaderInfo header = read_bam_header(bam.string());

  BuildOptions exclude;
  exclude.ambiguous_read_bases = AmbiguousBasePolicy::Exclude;
  RunIndex rejected;
  CHECK_THROWS_AS(
      covidx::detail::build_single_unit(header, exclude, 0, rejected),
      covidx::InputError);

  tests::TempDir temp;
  std::vector<tests::ReadSpec> reads = tests::tiny_reads();
  std::erase_if(
      reads, [](const tests::ReadSpec& read) { return read.name == "r_star"; });
  const fs::path no_star = temp.path() / "no_star.bam";
  tests::make_bam(no_star, tests::kTinyHeader, std::move(reads));
  const BamHeaderInfo no_star_header = read_bam_header(no_star.string());
  const std::vector<bool> naive = naive_coverage(no_star_header, exclude, 0);
  RunIndex index;
  covidx::detail::build_single_unit(no_star_header, exclude, 0, index);
  compare_index(index, naive, 7);
}

TEST_CASE("single-unit build handles a contig with no reads") {
  const fs::path bam = fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam";
  const BamHeaderInfo header = read_bam_header(bam.string());
  RunIndex index;
  covidx::detail::build_single_unit(header, BuildOptions{}, 1, index);
  CHECK(index.runs() == 0);
  CHECK(index.total_covered() == 0);
  CHECK(index.length() == 500);
}

TEST_CASE("effective exclude mask includes supplementary by default") {
  BuildOptions options;
  CHECK(effective_exclude_flags(options) == 0x0F04);
  options.include_supplementary = true;
  CHECK(effective_exclude_flags(options) == 0x0704);
  options.exclude_flags = 0;
  CHECK(effective_exclude_flags(options) == 0);
  options.include_supplementary = false;
  CHECK(effective_exclude_flags(options) == 0x0800);
}
