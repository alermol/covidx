#pragma once
//
// Maximal covered interval [start, end), 0-based half-open.
//
#include <cstdint>

namespace covidx::detail {

struct Run {
  uint64_t start = 0;
  uint64_t end = 0;
};

[[nodiscard]] inline bool operator==(const Run& a, const Run& b) {
  return a.start == b.start && a.end == b.end;
}

} // namespace covidx::detail
