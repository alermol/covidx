#include "covidx/detail/text.h"

#include <cstddef>

#include "covidx/error.h"

namespace covidx::detail {

bool is_valid_utf8(std::string_view text) {
  const auto* bytes = reinterpret_cast<const unsigned char*>(text.data());
  std::size_t i = 0;
  while (i < text.size()) {
    const unsigned char c = bytes[i];
    if (c < 0x80) {
      ++i;
      continue;
    }
    std::size_t extra = 0;
    unsigned cp = 0;
    if ((c & 0xE0) == 0xC0) {
      extra = 1;
      cp = c & 0x1F;
    } else if ((c & 0xF0) == 0xE0) {
      extra = 2;
      cp = c & 0x0F;
    } else if ((c & 0xF8) == 0xF0) {
      extra = 3;
      cp = c & 0x07;
    } else {
      return false;
    }
    if (i + extra >= text.size())
      return false;
    for (std::size_t k = 1; k <= extra; ++k) {
      const unsigned char next = bytes[i + k];
      if ((next & 0xC0) != 0x80)
        return false;
      cp = (cp << 6) | (next & 0x3F);
    }
    if (extra == 1 && cp < 0x80)
      return false;
    if (extra == 2 && cp < 0x800)
      return false;
    if (extra == 3 && cp < 0x10000)
      return false;
    if (cp >= 0xD800 && cp <= 0xDFFF)
      return false;
    if (cp > 0x10FFFF)
      return false;
    i += extra + 1;
  }
  return true;
}

void validate_unit_name(std::string_view name, std::string_view kind,
                        bool allow_pooled) {
  if (name.empty())
    throw UsageError(std::string(kind) + " name must not be empty");
  if (name.size() > 255) {
    throw UsageError(std::string(kind) + " name is longer than 255 bytes: '" +
                     std::string(name) + "'");
  }
  if (!is_valid_utf8(name)) {
    throw UsageError(std::string(kind) + " name is not valid UTF-8");
  }
  if (name.find('\0') != std::string_view::npos ||
      name.find('\t') != std::string_view::npos ||
      name.find('\r') != std::string_view::npos ||
      name.find('\n') != std::string_view::npos) {
    throw UsageError(std::string(kind) + " name contains a control character");
  }
  if (!allow_pooled && name == "pooled") {
    throw UsageError(std::string(kind) + " name 'pooled' is reserved");
  }
}

bool valid_index_name(std::string_view name) {
  if (name.empty() || name.size() > 128 || name == "." || name == "..")
    return false;
  for (const char c : name) {
    const bool ok = (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
                    (c >= '0' && c <= '9') || c == '.' || c == '_' || c == '-';
    if (!ok)
      return false;
  }
  return true;
}

} // namespace covidx::detail
