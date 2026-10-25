#pragma once
//
// SHA-256 (FIPS 180-4), adapted from gcrank (MIT).
//
#include <cstddef>
#include <cstdint>
#include <string>

namespace covidx::detail {

class SHA256 {
public:
  SHA256() {
    reset();
  }

  void reset();
  void update(const void* data, std::size_t n);
  void finish(uint8_t out[32]);

private:
  void transform(const uint8_t block[64]);

  uint32_t h_[8];
  uint8_t buf_[64];
  std::size_t buf_len_ = 0;
  uint64_t total_bits_ = 0;
};

// Lowercase hex digest of a memory buffer.
[[nodiscard]] std::string sha256_hex(const void* data, std::size_t n);

// Lowercase hex digest of the exact bytes of a file. Throws IoError.
[[nodiscard]] std::string sha256_file_hex(const std::string& path);

} // namespace covidx::detail
