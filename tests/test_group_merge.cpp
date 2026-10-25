#include "doctest.h"

#include <algorithm>
#include <cstdint>
#include <random>
#include <vector>

#include "covidx/detail/merge.h"

using covidx::detail::merge_run_lists;
using covidx::detail::Run;

namespace {

std::vector<Run> naive_union(const std::vector<std::vector<Run>>& lists) {
  std::vector<Run> runs;
  for (const auto& list : lists)
    runs.insert(runs.end(), list.begin(), list.end());
  std::sort(runs.begin(), runs.end(),
            [](const Run& a, const Run& b) { return a.start < b.start; });
  std::vector<Run> out;
  for (const Run& run : runs) {
    if (out.empty() || run.start > out.back().end) {
      out.push_back(run);
    } else if (run.end > out.back().end) {
      out.back().end = run.end;
    }
  }
  return out;
}

void check_equal(const std::vector<Run>& got, const std::vector<Run>& want) {
  REQUIRE(got.size() == want.size());
  for (std::size_t i = 0; i < want.size(); ++i) {
    CHECK(got[i].start == want[i].start);
    CHECK(got[i].end == want[i].end);
  }
}

std::vector<Run> random_canonical(std::mt19937& rng, uint64_t length) {
  std::vector<Run> runs;
  uint64_t cursor = 0;
  while (cursor < length) {
    cursor += rng() % 10;
    if (cursor >= length)
      break;
    const uint64_t end = std::min<uint64_t>(length, cursor + 1 + rng() % 20);
    runs.push_back(Run{cursor, end});
    cursor = end + 1 + rng() % 5;
  }
  return runs;
}

} // namespace

TEST_CASE("merge_run_lists merges overlaps and touching runs") {
  const std::vector<Run> a = {{0, 5}, {10, 20}, {30, 40}};
  const std::vector<Run> b = {{3, 12}, {20, 25}, {50, 60}};
  const std::vector<Run> c = {{18, 19}, {40, 50}};
  const std::vector<Run> got = merge_run_lists({&a, &b, &c});
  check_equal(got, naive_union({a, b, c}));
  REQUIRE(got.size() == 2);
  CHECK(got[0].start == 0);
  CHECK(got[0].end == 25);
  CHECK(got[1].start == 30);
  CHECK(got[1].end == 60);
}

TEST_CASE("merge_run_lists handles empty inputs") {
  const std::vector<Run> empty;
  const std::vector<Run> one = {{5, 10}};
  CHECK(merge_run_lists({}).empty());
  CHECK(merge_run_lists({&empty}).empty());
  check_equal(merge_run_lists({&empty, &one, &empty}), one);
}

TEST_CASE("merge_run_lists is order-independent") {
  const std::vector<Run> a = {{0, 5}, {20, 30}};
  const std::vector<Run> b = {{3, 25}, {40, 50}};
  check_equal(merge_run_lists({&a, &b}), merge_run_lists({&b, &a}));
}

TEST_CASE("merge_run_lists property test") {
  std::mt19937 rng(424242);
  for (int trial = 0; trial < 300; ++trial) {
    const std::size_t count = 1 + rng() % 6;
    std::vector<std::vector<Run>> lists;
    lists.reserve(count);
    for (std::size_t i = 0; i < count; ++i)
      lists.push_back(random_canonical(rng, 500));
    std::vector<const std::vector<Run>*> pointers;
    for (const auto& list : lists)
      pointers.push_back(&list);
    check_equal(merge_run_lists(pointers), naive_union(lists));
  }
}
