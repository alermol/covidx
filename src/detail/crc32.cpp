#include "covidx/detail/crc32.h"

#include <algorithm>
#include <fstream>
#include <vector>

#include <zlib.h>

#include "covidx/error.h"

namespace covidx::detail {

namespace {
constexpr std::size_t kBufferBytes = 1u << 18;
} // namespace

void CRC32::update(const void* data, std::size_t n) {
#if ZLIB_VERNUM >= 0x1290
  value_ = static_cast<uint32_t>(::crc32_z(
               static_cast<uLong>(value_ ^ 0xFFFFFFFFu),
               static_cast<const Bytef*>(data), static_cast<z_size_t>(n))) ^
           0xFFFFFFFFu;
#else
  uLong crc = static_cast<uLong>(value_ ^ 0xFFFFFFFFu);
  const auto* p = static_cast<const uint8_t*>(data);
  while (n > 0) {
    const uInt chunk = static_cast<uInt>(std::min<std::size_t>(n, 0xFFFFFFFFu));
    crc = ::crc32(crc, p, chunk);
    p += chunk;
    n -= chunk;
  }
  value_ = static_cast<uint32_t>(crc) ^ 0xFFFFFFFFu;
#endif
}

uint32_t crc32_file(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in)
    throw IoError("cannot open '" + path + "' for reading");
  CRC32 crc;
  std::vector<char> buf(kBufferBytes);
  while (in) {
    in.read(buf.data(), static_cast<std::streamsize>(buf.size()));
    const std::streamsize got = in.gcount();
    if (got > 0)
      crc.update(buf.data(), static_cast<std::size_t>(got));
  }
  if (in.bad())
    throw IoError("read error on '" + path + "'");
  return crc.value();
}

std::string crc32_hex(uint32_t crc) {
  static constexpr char kHex[] = "0123456789abcdef";
  std::string out(8, '0');
  for (int i = 7; i >= 0; --i) {
    out[static_cast<std::size_t>(i)] = kHex[crc & 0x0Fu];
    crc >>= 4;
  }
  return out;
}

} // namespace covidx::detail
