#include "covidx/detail/run_index.h"

#include <string>

#include "covidx/error.h"
#include "covidx/format.h"

namespace covidx::detail {

void RunIndex::reset(uint64_t length) {
  if (length > kMaxContigLength) {
    throw UsageError("contig length " + std::to_string(length) +
                     " exceeds BAM alignment coordinate limit INT32_MAX");
  }
  length_ = length;
  w_ = 20;
  total_covered_ = 0;
  runs_.clear();
  start_.clear();
  end_.clear();
  cum_.assign(1, 0);
  first_.clear();
}

void RunIndex::set_runs(std::vector<Run> runs) {
  uint64_t total = 0;
  for (std::size_t i = 0; i < runs.size(); ++i) {
    const Run& run = runs[i];
    if (run.start >= run.end) {
      throw FormatError("run " + std::to_string(i) + " is empty or reversed");
    }
    if (run.end > length_) {
      throw FormatError("run " + std::to_string(i) +
                        " ends beyond contig length " +
                        std::to_string(length_));
    }
    if (i > 0 && runs[i - 1].end >= run.start) {
      throw FormatError("runs " + std::to_string(i - 1) + " and " +
                        std::to_string(i) + " overlap or touch");
    }
    total += run.end - run.start;
  }
  runs_ = std::move(runs);
  start_.resize(runs_.size());
  end_.resize(runs_.size());
  cum_.assign(runs_.size() + 1, 0);
  for (std::size_t i = 0; i < runs_.size(); ++i) {
    start_[i] = static_cast<uint32_t>(runs_[i].start);
    end_[i] = static_cast<uint32_t>(runs_[i].end);
    cum_[i + 1] = cum_[i] + (runs_[i].end - runs_[i].start);
  }
  total_covered_ = total;
  w_ = choose_w(length_, runs_.size());
  first_ = build_first(runs_, length_, w_);
  validate();
}

void RunIndex::validate() const {
  if (length_ > kMaxContigLength)
    throw FormatError("contig length exceeds INT32_MAX");
  if (start_.size() != runs_.size() || end_.size() != runs_.size()) {
    throw FormatError("run array sizes disagree");
  }
  if (cum_.size() != runs_.size() + 1)
    throw FormatError("cum array size is not R+1");
  if (w_ < 6 || w_ > 20)
    throw FormatError("bucket width exponent is out of range");
  if (runs_.empty()) {
    if (!first_.empty())
      throw FormatError("empty run set has a bucket table");
    if (cum_.size() != 1 || cum_[0] != 0)
      throw FormatError("empty run set has a bad cum array");
    if (total_covered_ != 0)
      throw FormatError("empty run set has non-zero coverage");
    return;
  }
  if (runs_.size() > (length_ + 1) / 2) {
    throw FormatError("run count exceeds the canonical maximum");
  }
  if (first_.size() != bucket_entries(length_, w_)) {
    throw FormatError("bucket table size disagrees with length and width");
  }

  uint64_t cumulative = 0;
  for (std::size_t i = 0; i < runs_.size(); ++i) {
    const uint64_t s = runs_[i].start;
    const uint64_t e = runs_[i].end;
    if (s >= e || e > length_)
      throw FormatError("run is out of range");
    if (i > 0 && runs_[i - 1].end >= s)
      throw FormatError("runs overlap or touch");
    cumulative += e - s;
    if (cum_[i] != cumulative - (e - s))
      throw FormatError("cum recurrence is broken");
  }
  if (cum_.back() != cumulative || total_covered_ != cumulative) {
    throw FormatError("cum sentinel or total coverage is wrong");
  }
  if (cum_.back() > length_)
    throw FormatError("total coverage exceeds contig length");
  if (start_.front() != static_cast<uint32_t>(runs_.front().start) ||
      start_.back() != static_cast<uint32_t>(runs_.back().start) ||
      end_.front() != static_cast<uint32_t>(runs_.front().end) ||
      end_.back() != static_cast<uint32_t>(runs_.back().end)) {
    throw FormatError("uint32 run arrays disagree with runs");
  }

  if (first_.front() != 0 || first_.back() != runs_.size()) {
    throw FormatError("bucket table endpoints are wrong");
  }
  uint32_t previous = 0;
  uint32_t expected = 0;
  for (std::size_t b = 0; b < first_.size(); ++b) {
    if (first_[b] < previous || first_[b] > runs_.size()) {
      throw FormatError("bucket table is not monotone");
    }
    previous = first_[b];
    const uint64_t limit = static_cast<uint64_t>(b) << w_;
    while (expected < runs_.size() && runs_[expected].start < limit)
      ++expected;
    if (first_[b] != expected)
      throw FormatError("bucket table value is wrong");
  }
}

RunIndexView RunIndex::view() const {
  RunIndexView out;
  out.length = length_;
  out.runs = runs_.size();
  out.total_covered = total_covered_;
  out.w = w_;
  out.start = start_;
  out.end = end_;
  out.cum = cum_;
  out.first = first_;
  return out;
}

uint64_t RunIndex::covered_count(uint64_t x, uint64_t y) const {
  return detail::covered_count(view(), x, y);
}

} // namespace covidx::detail
