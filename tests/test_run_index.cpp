#include "doctest.h"

#include <algorithm>
#include <cstdint>
#include <random>
#include <vector>

#include "covidx/detail/run_index.h"
#include "covidx/error.h"
#include "covidx/format.h"

using covidx::detail::Run;
using covidx::detail::RunIndex;

TEST_CASE("RunIndex empty and small contigs") {
  RunIndex index;
  index.reset(0);
  index.set_runs({});
  CHECK(index.runs() == 0);
  CHECK(index.total_covered() == 0);
  CHECK(index.length() == 0);
  CHECK_THROWS_AS((void)index.covered_count(0, 1), covidx::QueryError);

  RunIndex one;
  one.reset(1);
  one.set_runs({{0, 1}});
  CHECK(one.runs() == 1);
  CHECK(one.total_covered() == 1);
  CHECK(one.covered_count(0, 1) == 1);
  CHECK(one.view().w == 6);

  RunIndex zero;
  zero.reset(0);
  CHECK_THROWS_AS(zero.set_runs({{0, 1}}), covidx::FormatError);
}

TEST_CASE("RunIndex builds cum and first arrays") {
  RunIndex index;
  index.reset(100);
  index.set_runs({{10, 20}, {30, 35}, {80, 100}});
  const auto view = index.view();
  REQUIRE(view.cum.size() == 4);
  CHECK(view.cum[0] == 0);
  CHECK(view.cum[1] == 10);
  CHECK(view.cum[2] == 15);
  CHECK(view.cum[3] == 35);
  CHECK(view.total_covered == 35);
  CHECK(view.runs == 3);
  CHECK(view.first.front() == 0);
  CHECK(view.first.back() == 3);
  index.validate();
}

TEST_CASE("RunIndex rejects non-canonical runs") {
  RunIndex index;
  index.reset(100);
  CHECK_THROWS_AS(index.set_runs({{20, 30}, {10, 20}}), covidx::FormatError);
  CHECK_THROWS_AS(index.set_runs({{10, 20}, {15, 25}}), covidx::FormatError);
  CHECK_THROWS_AS(index.set_runs({{10, 20}, {20, 30}}), covidx::FormatError);
  CHECK_THROWS_AS(index.set_runs({{10, 10}}), covidx::FormatError);
  CHECK_THROWS_AS(index.set_runs({{10, 9}}), covidx::FormatError);
  CHECK_THROWS_AS(index.set_runs({{90, 101}}), covidx::FormatError);
}

TEST_CASE("RunIndex length limit") {
  RunIndex index;
  index.reset(covidx::kMaxContigLength);
  index.set_runs({{covidx::kMaxContigLength - 1, covidx::kMaxContigLength}});
  CHECK(index.total_covered() == 1);
  CHECK_THROWS_AS(index.reset(covidx::kMaxContigLength + 1),
                  covidx::UsageError);
}

TEST_CASE("RunIndex query matches a naive bitmap") {
  std::mt19937 rng(998877);
  for (int trial = 0; trial < 30; ++trial) {
    const uint64_t length = 1 + rng() % 300;
    std::vector<bool> bitmap(length, false);
    for (uint64_t i = 0; i < length; ++i)
      bitmap[i] = (rng() % 4) == 0;
    std::vector<Run> runs;
    uint64_t i = 0;
    while (i < length) {
      if (!bitmap[i]) {
        ++i;
        continue;
      }
      const uint64_t start = i;
      while (i < length && bitmap[i])
        ++i;
      runs.push_back(Run{start, i});
    }
    RunIndex index;
    index.reset(length);
    index.set_runs(runs);
    CHECK(index.total_covered() == static_cast<uint64_t>(std::count(
                                       bitmap.begin(), bitmap.end(), true)));
    for (int q = 0; q < 300; ++q) {
      const uint64_t x = rng() % length;
      const uint64_t y = x + 1 + rng() % (length - x);
      uint64_t want = 0;
      for (uint64_t p = x; p < y; ++p)
        want += bitmap[p] ? 1 : 0;
      CHECK(index.covered_count(x, y) == want);
    }
  }
}

TEST_CASE("RunIndex dense alternating runs") {
  RunIndex index;
  index.reset(1001);
  std::vector<Run> runs;
  for (uint64_t i = 0; i < 1000; i += 2)
    runs.push_back(Run{i, i + 1});
  index.set_runs(runs);
  CHECK(index.runs() == 500);
  CHECK(index.total_covered() == 500);
  CHECK(index.covered_count(0, 1000) == 500);
  CHECK(index.covered_count(0, 1001) == 500);
  CHECK(index.covered_count(1, 2) == 0);
  index.validate();
}
