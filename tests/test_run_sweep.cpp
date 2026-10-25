#include "doctest.h"

#include <algorithm>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

#include "covidx/detail/run_sweep.h"

using covidx::detail::Block;
using covidx::detail::BlockBuffer;
using covidx::detail::Run;
using covidx::detail::RunSweep;

namespace {

std::vector<Run> canonical_union(std::vector<Block> blocks) {
  std::sort(blocks.begin(), blocks.end(), [](const Block& a, const Block& b) {
    return a.start < b.start || (a.start == b.start && a.end < b.end);
  });
  std::vector<Run> out;
  for (const Block& block : blocks) {
    if (out.empty() || block.start > out.back().end) {
      out.push_back(Run{block.start, block.end});
    } else if (block.end > out.back().end) {
      out.back().end = block.end;
    }
  }
  return out;
}

} // namespace

TEST_CASE("RunSweep merges touching and overlapping blocks") {
  std::vector<Run> runs;
  RunSweep sweep([&runs](Run run) { runs.push_back(run); });
  sweep.add(Block{10, 20});
  sweep.add(Block{0, 5});
  sweep.add(Block{5, 10});
  sweep.add(Block{30, 40});
  sweep.finish();
  REQUIRE(runs.size() == 2);
  CHECK(runs[0].start == 0);
  CHECK(runs[0].end == 20);
  CHECK(runs[1].start == 30);
  CHECK(runs[1].end == 40);
}

TEST_CASE("RunSweep respects a strict watermark") {
  std::vector<Run> runs;
  RunSweep sweep([&runs](Run run) { runs.push_back(run); });
  sweep.add(Block{10, 20});
  sweep.add(Block{0, 5});
  CHECK(runs.empty());
  sweep.advance(10);
  // The active [0,5) run is only emitted when the next non-adjacent block is
  // processed; it cannot be emitted before the watermark passes that block.
  CHECK(runs.empty());
  sweep.advance(20);
  REQUIRE(runs.size() == 1);
  CHECK(runs[0].start == 0);
  CHECK(runs[0].end == 5);
  sweep.finish();
  REQUIRE(runs.size() == 2);
  CHECK(runs[1].start == 10);
  CHECK(runs[1].end == 20);
}

TEST_CASE("RunSweep handles same-POS and nested blocks") {
  std::vector<Run> runs;
  RunSweep sweep([&runs](Run run) { runs.push_back(run); });
  sweep.add(Block{100, 110});
  sweep.add(Block{100, 105});
  sweep.add(Block{100, 120});
  sweep.finish();
  REQUIRE(runs.size() == 1);
  CHECK(runs[0].start == 100);
  CHECK(runs[0].end == 120);
}

TEST_CASE("RunSweep creates separate runs for one-base gaps") {
  std::vector<Run> runs;
  RunSweep sweep([&runs](Run run) { runs.push_back(run); });
  sweep.add(Block{0, 5});
  sweep.add(Block{6, 10});
  sweep.finish();
  REQUIRE(runs.size() == 2);
}

TEST_CASE("BlockBuffer drains in start order") {
  BlockBuffer buffer;
  buffer.add(Block{10, 20});
  buffer.add(Block{0, 5});
  buffer.add(Block{5, 8});
  std::vector<Block> drained;
  buffer.drain_before(10, [&drained](Block b) { drained.push_back(b); });
  REQUIRE(drained.size() == 2);
  CHECK(drained[0].start == 0);
  CHECK(drained[1].start == 5);
  CHECK(buffer.size() == 1);
  buffer.drain_all([&drained](Block b) { drained.push_back(b); });
  REQUIRE(drained.size() == 3);
  CHECK(drained[2].start == 10);
  CHECK(buffer.empty());
}

TEST_CASE("RunSweep matches a naive union for random blocks") {
  std::mt19937 rng(1234567);
  std::uniform_int_distribution<uint64_t> start_dist(0, 5000);
  std::uniform_int_distribution<uint64_t> len_dist(1, 50);

  for (int trial = 0; trial < 200; ++trial) {
    std::vector<Block> blocks;
    const int count = static_cast<int>(rng() % 40);
    for (int i = 0; i < count; ++i) {
      const uint64_t start = start_dist(rng);
      blocks.push_back(Block{start, start + len_dist(rng)});
    }
    std::vector<Run> runs;
    RunSweep sweep([&runs](Run run) { runs.push_back(run); });
    for (const Block& block : blocks)
      sweep.add(block);
    sweep.finish();
    const std::vector<Run> want = canonical_union(blocks);
    REQUIRE(runs.size() == want.size());
    for (std::size_t i = 0; i < want.size(); ++i) {
      CHECK(runs[i].start == want[i].start);
      CHECK(runs[i].end == want[i].end);
    }
  }
}
