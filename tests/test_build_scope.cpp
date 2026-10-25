#include "doctest.h"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <random>
#include <utility>
#include <vector>

#include "covidx/cov_index.h"
#include "covidx/detail/bam_reader.h"
#include "covidx/detail/build.h"
#include "covidx/detail/run_index.h"
#include "covidx/detail/scope.h"
#include "covidx/error.h"
#include "covidx/manifest.h"
#include "test_util.h"

namespace fs = std::filesystem;
using covidx::BuildOptions;
using covidx::Region;
using covidx::detail::BamHeaderInfo;
using covidx::detail::Block;
using covidx::detail::clip_blocks;
using covidx::detail::make_scope;
using covidx::detail::read_bam_header;
using covidx::detail::RunIndex;
using covidx::detail::Scope;
using covidx::detail::window_in_scope;

namespace {

BamHeaderInfo tiny_header() {
  return read_bam_header(
      (fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam").string());
}

bool runs_in_scope(const RunIndex& index,
                   const std::vector<std::pair<uint64_t, uint64_t>>& scope) {
  for (const auto& run : index.run_list()) {
    bool inside = false;
    for (const auto& [start, end] : scope) {
      if (run.start >= start && run.end <= end) {
        inside = true;
        break;
      }
    }
    if (!inside)
      return false;
  }
  return true;
}

} // namespace

TEST_CASE("default scope covers every contig") {
  const BamHeaderInfo header = tiny_header();
  const Scope scope = make_scope(header, BuildOptions{});
  REQUIRE(scope.intervals.size() == 2);
  CHECK(scope.intervals[0].contig == "chr1");
  CHECK(scope.intervals[0].start == 0);
  CHECK(scope.intervals[0].end == 1000);
  CHECK(scope.intervals[1].contig == "chr2");
  CHECK(scope.intervals[1].end == 500);
  CHECK_FALSE(scope.partial);
}

TEST_CASE("--contigs selects full contigs and deduplicates") {
  const BamHeaderInfo header = tiny_header();
  BuildOptions options;
  options.contigs = {"chr2", "chr2"};
  const Scope scope = make_scope(header, options);
  REQUIRE(scope.intervals.size() == 1);
  CHECK(scope.intervals[0].contig == "chr2");
  CHECK(scope.intervals[0].start == 0);
  CHECK(scope.intervals[0].end == 500);
  CHECK_FALSE(scope.partial);

  options.contigs = {"nope"};
  CHECK_THROWS_AS((void)make_scope(header, options), covidx::UsageError);
}

TEST_CASE("--regions sorts, merges, and marks partial scope") {
  const BamHeaderInfo header = tiny_header();
  BuildOptions options;
  options.regions = {
      Region{"chr1", 340, 650},
      Region{"chr1", 100, 350},
      Region{"chr1", 900, 1000},
      Region{"chr2", 0, 500},
  };
  const Scope scope = make_scope(header, options);
  REQUIRE(scope.intervals.size() == 3);
  CHECK(scope.intervals[0].contig == "chr1");
  CHECK(scope.intervals[0].start == 100);
  CHECK(scope.intervals[0].end == 650);
  CHECK(scope.intervals[1].start == 900);
  CHECK(scope.intervals[1].end == 1000);
  CHECK(scope.intervals[2].contig == "chr2");
  CHECK(scope.intervals[2].start == 0);
  CHECK(scope.intervals[2].end == 500);
  CHECK(scope.partial);
}

TEST_CASE("full-length --regions are not partial") {
  const BamHeaderInfo header = tiny_header();
  BuildOptions options;
  options.regions = {Region{"chr1", 0, 1000}, Region{"chr2", 0, 500}};
  const Scope scope = make_scope(header, options);
  CHECK_FALSE(scope.partial);
  CHECK(scope.intervals.size() == 2);
}

TEST_CASE("scope validation errors") {
  const BamHeaderInfo header = tiny_header();
  BuildOptions options;
  options.contigs = {"chr1"};
  options.regions = {Region{"chr1", 0, 10}};
  CHECK_THROWS_AS((void)make_scope(header, options), covidx::UsageError);

  BuildOptions unknown;
  unknown.regions = {Region{"nope", 0, 10}};
  CHECK_THROWS_AS((void)make_scope(header, unknown), covidx::InputError);

  BuildOptions empty;
  empty.regions = {Region{"chr1", 10, 10}};
  CHECK_THROWS_AS((void)make_scope(header, empty), covidx::InputError);

  BuildOptions beyond;
  beyond.regions = {Region{"chr1", 0, 1001}};
  CHECK_THROWS_AS((void)make_scope(header, beyond), covidx::InputError);
}

TEST_CASE("clip_blocks splits at scope gaps and discards outside blocks") {
  std::vector<std::pair<uint64_t, uint64_t>> intervals = {
      {0, 3}, {5, 10}, {20, 25}};
  const std::vector<Block> clipped = clip_blocks(
      {Block{0, 10}, Block{12, 15}, Block{19, 22}, Block{24, 30}}, intervals);
  REQUIRE(clipped.size() == 4);
  CHECK(clipped[0].start == 0);
  CHECK(clipped[0].end == 3);
  CHECK(clipped[1].start == 5);
  CHECK(clipped[1].end == 10);
  CHECK(clipped[2].start == 20);
  CHECK(clipped[2].end == 22);
  CHECK(clipped[3].start == 24);
  CHECK(clipped[3].end == 25);

  const std::vector<Block> merged =
      clip_blocks({Block{0, 10}}, {{0, 5}, {5, 10}});
  REQUIRE(merged.size() == 1);
  CHECK(merged[0].start == 0);
  CHECK(merged[0].end == 10);
}

TEST_CASE("window_in_scope requires containment in one interval") {
  const std::vector<Region> intervals = {Region{"chr1", 100, 200},
                                         Region{"chr1", 300, 400}};
  CHECK(window_in_scope(intervals, 100, 200));
  CHECK(window_in_scope(intervals, 120, 180));
  CHECK_FALSE(window_in_scope(intervals, 90, 120));
  CHECK_FALSE(window_in_scope(intervals, 180, 220));
  CHECK_FALSE(window_in_scope(intervals, 200, 300));
  CHECK_FALSE(window_in_scope(intervals, 399, 401));
  CHECK_FALSE(window_in_scope(intervals, 150, 150));
  CHECK_FALSE(window_in_scope({}, 0, 1));
}

TEST_CASE("scoped build equals a full build restricted to the scope") {
  const BamHeaderInfo header = tiny_header();
  RunIndex full;
  covidx::detail::build_single_unit(header, BuildOptions{}, 0, full);

  const std::vector<std::pair<uint64_t, uint64_t>> scope = {
      {100, 350}, {500, 650}, {900, 1000}};
  RunIndex scoped;
  covidx::detail::build_single_unit_scoped(header, BuildOptions{}, 0, scope,
                                           scoped);
  CHECK(runs_in_scope(scoped, scope));

  std::mt19937 rng(31415);
  for (int trial = 0; trial < 500; ++trial) {
    const auto& interval = scope[rng() % scope.size()];
    const uint64_t x =
        interval.first + rng() % (interval.second - interval.first);
    const uint64_t y = x + 1 + rng() % (interval.second - x);
    CHECK(scoped.covered_count(x, y) == full.covered_count(x, y));
  }
  for (const auto& [start, end] : scope) {
    CHECK(scoped.covered_count(start, end) == full.covered_count(start, end));
  }
}

TEST_CASE("a read clipped away by scope is counted as no_blocks") {
  tests::TempDir temp;
  const fs::path bam = fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam";
  BuildOptions options;
  // r_d is 5M2D5M at chr1:200, so [205,206) lies entirely inside its
  // uncounted deletion. The read is fetched because its span overlaps the
  // scope, but clipping leaves no covered block.
  options.regions = {Region{"chr1", 205, 206}};
  options.count_deletions = false;
  options.work_dir = temp.path().string();
  covidx::CovIndex index;
  index.build(bam.string(), options);
  CHECK(index.build_stats().records_examined == 1);
  CHECK(index.build_stats().records_no_blocks == 1);
  CHECK(index.build_stats().records_accepted == 0);
  CHECK(index.manifest().source.records_accepted == 0);
  CHECK(index.covered_count("S2", "chr1", 205, 206) == 0);
}

TEST_CASE("random scopes match a full build restricted to the same intervals") {
  const BamHeaderInfo header = tiny_header();
  RunIndex full;
  covidx::detail::build_single_unit(header, BuildOptions{}, 0, full);
  std::mt19937 rng(271828);
  for (int trial = 0; trial < 20; ++trial) {
    std::vector<std::pair<uint64_t, uint64_t>> scope;
    uint64_t cursor = 0;
    while (cursor < 950) {
      cursor += 1 + rng() % 60;
      const uint64_t end = std::min<uint64_t>(1000, cursor + 1 + rng() % 120);
      if (cursor < end)
        scope.push_back({cursor, end});
      cursor = end + 1 + rng() % 40;
    }
    RunIndex scoped;
    covidx::detail::build_single_unit_scoped(header, BuildOptions{}, 0, scope,
                                             scoped);
    CHECK(runs_in_scope(scoped, scope));
    for (int q = 0; q < 100; ++q) {
      const auto& interval = scope[rng() % scope.size()];
      const uint64_t x =
          interval.first + rng() % (interval.second - interval.first);
      const uint64_t y = x + 1 + rng() % (interval.second - x);
      CHECK(scoped.covered_count(x, y) == full.covered_count(x, y));
    }
  }
}
