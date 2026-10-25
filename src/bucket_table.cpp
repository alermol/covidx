#include "covidx/detail/bucket_table.h"

#include <algorithm>
#include <string>

#include "covidx/error.h"
#include "covidx/format.h"

namespace covidx::detail {

uint32_t choose_w(uint64_t length, uint64_t runs) {
  if (length > kMaxContigLength) {
    throw UsageError("contig length " + std::to_string(length) +
                     " exceeds BAM alignment coordinate limit INT32_MAX");
  }
  if (runs == 0)
    return 20;
  uint32_t w = 6;
  while (w < 20 && (length >> w) > runs)
    ++w;
  return w;
}

uint64_t bucket_entries(uint64_t length, uint32_t w) {
  if (w > 20)
    throw UsageError("bucket width exponent exceeds 20");
  const uint64_t width = uint64_t{1} << w;
  return (length + width - 1) / width + 1;
}

std::vector<uint32_t> build_first(const std::vector<Run>& runs, uint64_t length,
                                  uint32_t w) {
  if (runs.empty())
    return {};
  const uint64_t entries = bucket_entries(length, w);
  std::vector<uint32_t> first(static_cast<std::size_t>(entries), 0);
  const uint64_t width = uint64_t{1} << w;
  std::size_t j = 0;
  for (uint64_t b = 0; b < entries; ++b) {
    const uint64_t limit = b * width;
    while (j < runs.size() && runs[j].start < limit)
      ++j;
    first[static_cast<std::size_t>(b)] = static_cast<uint32_t>(j);
  }
  return first;
}

uint64_t prefix_covered(const RunIndexView& view, uint64_t x) {
  if (x == 0 || view.runs == 0)
    return 0;
  if (x >= view.length)
    return view.total_covered;
  if (view.first.size() < 2)
    throw FormatError("bucket table is too short");

  const uint64_t b = x >> view.w;
  if (b + 1 >= view.first.size())
    throw FormatError("bucket index is out of range");
  const uint64_t lo = view.first[static_cast<std::size_t>(b)];
  const uint64_t hi = view.first[static_cast<std::size_t>(b + 1)];
  if (lo > hi || hi > view.start.size())
    throw FormatError("bucket table is not monotone");

  const auto begin = view.start.begin() + static_cast<std::ptrdiff_t>(lo);
  const auto end = view.start.begin() + static_cast<std::ptrdiff_t>(hi);
  const auto it = std::upper_bound(begin, end, x);
  const uint64_t n = static_cast<uint64_t>(it - view.start.begin());
  if (n == 0)
    return 0;

  const uint64_t j = n - 1;
  const uint64_t run_start = view.start[static_cast<std::size_t>(j)];
  const uint64_t run_end = view.end[static_cast<std::size_t>(j)];
  return view.cum[static_cast<std::size_t>(j)] + std::min(x, run_end) -
         run_start;
}

uint64_t covered_count(const RunIndexView& view, uint64_t x, uint64_t y) {
  if (x >= y)
    throw QueryError("empty or reversed window [" + std::to_string(x) + ", " +
                     std::to_string(y) + ")");
  if (y > view.length) {
    throw QueryError("window end " + std::to_string(y) +
                     " exceeds contig length " + std::to_string(view.length));
  }
  return prefix_covered(view, y) - prefix_covered(view, x);
}

} // namespace covidx::detail
