#pragma once
//
// Coarse bucket table over canonical run starts and the O(1) prefix-count.
//
#include <cstdint>
#include <span>
#include <vector>

#include "covidx/detail/run.h"

namespace covidx::detail {

// Smallest w in [6,20] with floor(length / 2^w) <= runs. Returns 20 when
// runs == 0. Throws UsageError when length > INT32_MAX.
[[nodiscard]] uint32_t choose_w(uint64_t length, uint64_t runs);

// Number of bucket entries B = ceil(length / 2^w) + 1.
[[nodiscard]] uint64_t bucket_entries(uint64_t length, uint32_t w);

// first[b] = number of runs with start < b * 2^w.
[[nodiscard]] std::vector<uint32_t> build_first(const std::vector<Run>& runs,
                                                uint64_t length, uint32_t w);

// Non-owning view of the arrays in a unit file. All positions are uint32;
// cum is uint64 with R+1 entries; first has B entries.
struct RunIndexView {
  uint64_t length = 0;
  uint64_t runs = 0;
  uint64_t total_covered = 0;
  uint32_t w = 20;
  std::span<const uint32_t> start;
  std::span<const uint32_t> end;
  std::span<const uint64_t> cum;
  std::span<const uint32_t> first;
};

// Number of covered positions in [0, x).
[[nodiscard]] uint64_t prefix_covered(const RunIndexView& view, uint64_t x);

// Number of covered positions in [x, y). Throws QueryError for x >= y or
// y > view.length.
[[nodiscard]] uint64_t covered_count(const RunIndexView& view, uint64_t x,
                                     uint64_t y);

} // namespace covidx::detail
