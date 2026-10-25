#include "covidx/detail/run_sweep.h"

#include <utility>

namespace covidx::detail {

void BlockBuffer::add(Block block) {
  pending_.push(block);
}

RunSweep::RunSweep(RunSink sink) : sink_(std::move(sink)) {}

void RunSweep::add(Block block) {
  buffer_.add(block);
}

void RunSweep::merge(const Block& block) {
  if (!active_) {
    active_ = true;
    cs_ = block.start;
    ce_ = block.end;
    return;
  }
  if (block.start <= ce_) {
    if (block.end > ce_)
      ce_ = block.end;
    return;
  }
  flush();
  active_ = true;
  cs_ = block.start;
  ce_ = block.end;
}

void RunSweep::flush() {
  if (!active_)
    return;
  if (ce_ > cs_)
    sink_(Run{cs_, ce_});
  active_ = false;
}

void RunSweep::advance(uint64_t watermark) {
  buffer_.drain_before(watermark, [this](const Block& block) { merge(block); });
}

void RunSweep::finish() {
  buffer_.drain_all([this](const Block& block) { merge(block); });
  flush();
}

} // namespace covidx::detail
