#include "doctest.h"

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <random>
#include <vector>

#include "covidx/detail/bucket_table.h"
#include "covidx/error.h"
#include "covidx/format.h"

using covidx::detail::bucket_entries;
using covidx::detail::build_first;
using covidx::detail::choose_w;
using covidx::detail::covered_count;
using covidx::detail::prefix_covered;
using covidx::detail::Run;
using covidx::detail::RunIndexView;

namespace {

struct OwnedIndex {
  std::vector<uint32_t> start;
  std::vector<uint32_t> end;
  std::vector<uint64_t> cum;
  std::vector<uint32_t> first;
  RunIndexView view;
};

bool canonical(const std::vector<Run>& runs, uint64_t length) {
  for (std::size_t i = 0; i < runs.size(); ++i) {
    if (runs[i].start >= runs[i].end || runs[i].end > length)
      return false;
    if (i > 0 && runs[i - 1].end >= runs[i].start)
      return false;
  }
  return true;
}

OwnedIndex make_index(const std::vector<Run>& runs, uint64_t length) {
  REQUIRE(canonical(runs, length));
  OwnedIndex out;
  const uint32_t w = choose_w(length, runs.size());
  out.start.reserve(runs.size());
  out.end.reserve(runs.size());
  out.cum.resize(runs.size() + 1, 0);
  for (std::size_t i = 0; i < runs.size(); ++i) {
    out.start.push_back(static_cast<uint32_t>(runs[i].start));
    out.end.push_back(static_cast<uint32_t>(runs[i].end));
    out.cum[i + 1] = out.cum[i] + (runs[i].end - runs[i].start);
  }
  out.first = build_first(runs, length, w);
  out.view.length = length;
  out.view.runs = runs.size();
  out.view.total_covered = out.cum.back();
  out.view.w = w;
  out.view.start = out.start;
  out.view.end = out.end;
  out.view.cum = out.cum;
  out.view.first = out.first;
  return out;
}

std::vector<Run> bitmap_runs(const std::vector<bool>& covered) {
  std::vector<Run> runs;
  std::size_t i = 0;
  while (i < covered.size()) {
    if (!covered[i]) {
      ++i;
      continue;
    }
    const std::size_t start = i;
    while (i < covered.size() && covered[i])
      ++i;
    runs.push_back(Run{start, i});
  }
  return runs;
}

} // namespace

TEST_CASE("choose_w follows the clamped smallest-width rule") {
  CHECK(choose_w(0, 0) == 20);
  CHECK(choose_w(1000, 0) == 20);
  CHECK(choose_w(1000, 500) == 6);
  CHECK(choose_w(1000, 14) == 7);
  CHECK(choose_w(1000, 15) == 6);
  CHECK(choose_w(10, 1) == 6);
  CHECK(choose_w(64, 1) == 6);
  CHECK(choose_w(65, 1) == 6);
  CHECK(choose_w(covidx::kMaxContigLength, 1) == 20);
  CHECK_THROWS_AS((void)choose_w(covidx::kMaxContigLength + 1, 1),
                  covidx::UsageError);
}

TEST_CASE("bucket_entries computes ceil(L/W)+1") {
  CHECK(bucket_entries(10, 6) == 2);
  CHECK(bucket_entries(64, 6) == 2);
  CHECK(bucket_entries(65, 6) == 3);
  CHECK(bucket_entries(1000, 6) == 17);
  CHECK(bucket_entries(1024, 6) == 17);
  CHECK(bucket_entries(1025, 6) == 18);
  CHECK(bucket_entries(10, 20) == 2);
}

TEST_CASE("build_first counts run starts below each bucket boundary") {
  const uint64_t length = 200;
  const uint32_t w = 6;
  const std::vector<Run> runs = {
      {0, 1}, {63, 64}, {64, 65}, {65, 66}, {128, 129}};
  const std::vector<uint32_t> first = build_first(runs, length, w);
  REQUIRE(first.size() == 5);
  CHECK(first[0] == 0);
  CHECK(first[1] == 2);
  CHECK(first[2] == 4);
  CHECK(first[3] == 5);
  CHECK(first[4] == 5);
  CHECK(build_first({}, length, w).empty());
}

TEST_CASE("prefix_covered edge cases") {
  const OwnedIndex index = make_index({{10, 20}, {30, 35}, {80, 100}}, 100);
  CHECK(prefix_covered(index.view, 0) == 0);
  CHECK(prefix_covered(index.view, 9) == 0);
  CHECK(prefix_covered(index.view, 10) == 0);
  CHECK(prefix_covered(index.view, 11) == 1);
  CHECK(prefix_covered(index.view, 19) == 9);
  CHECK(prefix_covered(index.view, 20) == 10);
  CHECK(prefix_covered(index.view, 21) == 10);
  CHECK(prefix_covered(index.view, 30) == 10);
  CHECK(prefix_covered(index.view, 31) == 11);
  CHECK(prefix_covered(index.view, 35) == 15);
  CHECK(prefix_covered(index.view, 36) == 15);
  CHECK(prefix_covered(index.view, 79) == 15);
  CHECK(prefix_covered(index.view, 80) == 15);
  CHECK(prefix_covered(index.view, 100) == 35);
  CHECK(prefix_covered(index.view, 200) == 35);
}

TEST_CASE("covered_count windows") {
  const OwnedIndex index = make_index({{10, 20}, {30, 35}, {80, 100}}, 100);
  CHECK(covered_count(index.view, 10, 20) == 10);
  CHECK(covered_count(index.view, 0, 100) == 35);
  CHECK(covered_count(index.view, 0, 10) == 0);
  CHECK(covered_count(index.view, 20, 30) == 0);
  CHECK(covered_count(index.view, 15, 25) == 5);
  CHECK_THROWS_AS((void)covered_count(index.view, 5, 5), covidx::QueryError);
  CHECK_THROWS_AS((void)covered_count(index.view, 6, 5), covidx::QueryError);
  CHECK_THROWS_AS((void)covered_count(index.view, 0, 101), covidx::QueryError);
}

TEST_CASE("predecessor in an earlier bucket and empty buckets") {
  const OwnedIndex index = make_index({{1, 2}, {450, 451}}, 500);
  CHECK(prefix_covered(index.view, 1) == 0);
  CHECK(prefix_covered(index.view, 2) == 1);
  CHECK(prefix_covered(index.view, 300) == 1);
  CHECK(prefix_covered(index.view, 450) == 1);
  CHECK(prefix_covered(index.view, 451) == 2);
}

TEST_CASE("dense alternating runs") {
  std::vector<Run> runs;
  for (uint64_t i = 0; i < 1000; i += 2)
    runs.push_back(Run{i, i + 1});
  const OwnedIndex index = make_index(runs, 1000);
  for (uint64_t x = 0; x <= 1000; ++x) {
    CHECK(prefix_covered(index.view, x) == (x + 1) / 2);
  }
  CHECK(covered_count(index.view, 0, 1000) == 500);
  CHECK(covered_count(index.view, 0, 1) == 1);
  CHECK(covered_count(index.view, 1, 2) == 0);
}

TEST_CASE("empty run set and zero-length contig") {
  const OwnedIndex empty = make_index({}, 0);
  CHECK(empty.view.runs == 0);
  CHECK(prefix_covered(empty.view, 0) == 0);
  CHECK(prefix_covered(empty.view, 10) == 0);
  CHECK_THROWS_AS((void)covered_count(empty.view, 0, 1), covidx::QueryError);
}

TEST_CASE("prefix_covered matches a naive bitmap for random runs") {
  std::mt19937 rng(20261024);
  for (int trial = 0; trial < 50; ++trial) {
    const uint64_t length = 1 + rng() % 400;
    std::vector<bool> bitmap(length, false);
    for (uint64_t i = 0; i < length; ++i)
      bitmap[i] = (rng() % 3) == 0;
    const std::vector<Run> runs = bitmap_runs(bitmap);
    const OwnedIndex index = make_index(runs, length);

    uint64_t naive = 0;
    for (uint64_t x = 0; x <= length; ++x) {
      CHECK(prefix_covered(index.view, x) == naive);
      if (x < length && bitmap[x])
        ++naive;
    }
    CHECK(index.view.total_covered == naive);
    for (int window = 0; window < 200; ++window) {
      const uint64_t x = rng() % length;
      const uint64_t y = x + 1 + rng() % (length - x);
      uint64_t want = 0;
      for (uint64_t i = x; i < y; ++i)
        want += bitmap[i] ? 1 : 0;
      CHECK(covered_count(index.view, x, y) == want);
    }
  }
}
