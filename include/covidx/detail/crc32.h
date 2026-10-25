#pragma once
//
// CRC-32/ISO-HDLC (IEEE 802.3, reflected polynomial 0xEDB88320).
//
// Adapted from gcrank (MIT), pinned at
// 704eee87135ff221ec3c4d60068c39c48bd0bfaa.
//
#include <cstddef>
#include <cstdint>
#include <string>

namespace covidx::detail {

// Streaming CRC-32. A default-constructed object is ready for updates.
class CRC32 {
public:
  CRC32() = default;

  void update(const void* data, std::size_t n);
  [[nodiscard]] uint32_t value() const {
    return value_ ^ 0xFFFFFFFFu;
  }
  void reset() {
    value_ = 0xFFFFFFFFu;
  }

  [[nodiscard]] static uint32_t compute(const void* data, std::size_t n) {
    CRC32 crc;
    crc.update(data, n);
    return crc.value();
  }

private:
  uint32_t value_ = 0xFFFFFFFFu;
};

// CRC-32 of the exact bytes of a file. Throws IoError on read failure.
[[nodiscard]] uint32_t crc32_file(const std::string& path);

// Formats a checksum as exactly eight lowercase hexadecimal digits.
[[nodiscard]] std::string crc32_hex(uint32_t crc);

} // namespace covidx::detail
