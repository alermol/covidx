#include <cstddef>
#include <cstdint>
#include <string_view>

#include "covidx/detail/json.h"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, std::size_t size) {
  try {
    (void)covidx::detail::parse_json(
        std::string_view(reinterpret_cast<const char*>(data), size), "fuzz");
  } catch (...) {
  }
  return 0;
}
