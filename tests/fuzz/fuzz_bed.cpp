#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string>

#include "covidx/detail/bed.h"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, std::size_t size) {
  try {
    std::istringstream in(
        std::string(reinterpret_cast<const char*>(data), size));
    (void)covidx::detail::read_bed_stream(in, "fuzz");
  } catch (...) {
  }
  return 0;
}
