#include <cstddef>
#include <cstdint>
#include <vector>

#include "covidx/detail/run_index.h"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, std::size_t size) {
  if (size < 4)
    return 0;
  const uint64_t length = 1 + ((data[0] | (data[1] << 8)) % 500);
  covidx::detail::RunIndex index;
  index.reset(length);
  std::vector<covidx::detail::Run> runs;
  uint64_t position = 0;
  for (std::size_t i = 2; i + 1 < size && position < length; i += 2) {
    position += data[i] % 7;
    const uint64_t end = position + 1 + (data[i + 1] % 9);
    if (position >= end || end > length)
      break;
    runs.push_back(covidx::detail::Run{position, end});
    position = end + 1;
  }
  try {
    index.set_runs(std::move(runs));
    for (std::size_t i = 2; i + 3 < size; i += 4) {
      const uint64_t x = data[i] % length;
      const uint64_t y =
          x + 1 + ((data[i + 1] << 8 | data[i + 2]) % (length - x));
      (void)index.covered_count(x, y);
    }
  } catch (...) {
  }
  return 0;
}
