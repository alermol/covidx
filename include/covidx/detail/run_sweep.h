#pragma once
//
// Interval buffering and canonical run construction.
//
#include <cstddef>
#include <cstdint>
#include <functional>
#include <queue>
#include <vector>

#include "covidx/detail/cigar.h"
#include "covidx/detail/run.h"

namespace covidx::detail {

// A min-heap of blocks used to reorder blocks within one read or one unit
// sweep. Blocks are drained in (start, end) order.
class BlockBuffer {
public:
  void add(Block block);

  // Calls sink for every block with start < watermark, in order.
  template <class F> void drain_before(uint64_t watermark, F&& sink) {
    while (!pending_.empty() && pending_.top().start < watermark) {
      sink(pending_.top());
      pending_.pop();
    }
  }

  template <class F> void drain_all(F&& sink) {
    while (!pending_.empty()) {
      sink(pending_.top());
      pending_.pop();
    }
  }

  [[nodiscard]] std::size_t size() const {
    return pending_.size();
  }
  [[nodiscard]] bool empty() const {
    return pending_.empty();
  }

private:
  struct MinStart {
    bool operator()(const Block& a, const Block& b) const {
      if (a.start != b.start)
        return a.start > b.start;
      return a.end > b.end;
    }
  };
  std::priority_queue<Block, std::vector<Block>, MinStart> pending_;
};

using RunSink = std::function<void(Run)>;

// Merges blocks into canonical maximal runs. add() accepts blocks in any
// order; advance(watermark) finalizes runs that can no longer grow, and
// finish() drains everything.
class RunSweep {
public:
  explicit RunSweep(RunSink sink);

  void advance(uint64_t watermark);
  void add(Block block);
  void finish();

private:
  void merge(const Block& block);
  void flush();

  BlockBuffer buffer_;
  RunSink sink_;
  bool active_ = false;
  uint64_t cs_ = 0;
  uint64_t ce_ = 0;
};

} // namespace covidx::detail
