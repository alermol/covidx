#include "doctest.h"

#include <cstdint>
#include <vector>

#include "covidx/detail/bucket_table.h"
#include "covidx/format.h"

using covidx::detail::bucket_entries;
using covidx::detail::build_first;
using covidx::detail::choose_w;
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

OwnedIndex make_index(const std::vector<Run>& runs, uint64_t length) {
  OwnedIndex out;
  const uint32_t w = choose_w(length, runs.size());
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

} // namespace

TEST_CASE("bucket boundaries at W-1, W, and W+1") {
  const uint64_t length = 200;
  const uint32_t w = 6;
  const std::vector<Run> runs = {
      {63, 64}, {64, 65}, {65, 66}, {127, 128}, {128, 129}};
  const std::vector<uint32_t> first = build_first(runs, length, w);
  REQUIRE(first.size() == bucket_entries(length, w));
  CHECK(first[0] == 0);
  CHECK(first[1] == 1);
  CHECK(first[2] == 4);
  CHECK(first[3] == 5);
  CHECK(first[4] == 5);

  const OwnedIndex index = make_index(runs, length);
  CHECK(prefix_covered(index.view, 63) == 0);
  CHECK(prefix_covered(index.view, 64) == 1);
  CHECK(prefix_covered(index.view, 65) == 2);
  CHECK(prefix_covered(index.view, 66) == 3);
  CHECK(prefix_covered(index.view, 127) == 3);
  CHECK(prefix_covered(index.view, 128) == 4);
  CHECK(prefix_covered(index.view, 129) == 5);
}

TEST_CASE("length exact multiple of the bucket width") {
  const OwnedIndex index = make_index({{100, 128}}, 128);
  CHECK(index.view.w == 7);
  CHECK(index.view.first.size() == 2);
  CHECK(index.view.first.back() == 1);
  CHECK(prefix_covered(index.view, 100) == 0);
  CHECK(prefix_covered(index.view, 127) == 27);
  CHECK(prefix_covered(index.view, 128) == 28);
  CHECK(prefix_covered(index.view, 1000) == 28);
}

TEST_CASE("contig near INT32_MAX uses uint32 positions") {
  const uint64_t length = covidx::kMaxContigLength;
  const OwnedIndex index = make_index({{length - 1, length}}, length);
  CHECK(index.view.start.front() == static_cast<uint32_t>(length - 1));
  CHECK(index.view.end.front() == static_cast<uint32_t>(length));
  CHECK(index.view.w == 20);
  CHECK(index.view.first.size() == bucket_entries(length, 20));
  CHECK(index.view.first.back() == 1);
  CHECK(prefix_covered(index.view, length - 1) == 0);
  CHECK(prefix_covered(index.view, length) == 1);
  CHECK(prefix_covered(index.view, length + 100) == 1);
}

TEST_CASE("single bucket when the contig is narrower than W") {
  const OwnedIndex index = make_index({{1, 3}, {5, 7}}, 10);
  CHECK(index.view.w == 6);
  CHECK(index.view.first.size() == 2);
  CHECK(index.view.first[0] == 0);
  CHECK(index.view.first[1] == 2);
  CHECK(prefix_covered(index.view, 0) == 0);
  CHECK(prefix_covered(index.view, 1) == 0);
  CHECK(prefix_covered(index.view, 3) == 2);
  CHECK(prefix_covered(index.view, 10) == 4);
}
