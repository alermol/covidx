#include "doctest.h"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

#include "covidx/detail/crc32.h"
#include "covidx/detail/sha256.h"
#include "covidx/error.h"

namespace fs = std::filesystem;

TEST_CASE("CRC32 known vector") {
  const std::string text = "123456789";
  CHECK(covidx::detail::CRC32::compute(text.data(), text.size()) ==
        0xCBF43926u);
  CHECK(covidx::detail::crc32_hex(0xCBF43926u) == "cbf43926");
  CHECK(covidx::detail::crc32_hex(0) == "00000000");
}

TEST_CASE("CRC32 streaming equals one-shot") {
  const std::string text = "the quick brown fox jumps over the lazy dog";
  covidx::detail::CRC32 streaming;
  streaming.update(text.data(), 10);
  streaming.update(text.data() + 10, text.size() - 10);
  CHECK(streaming.value() ==
        covidx::detail::CRC32::compute(text.data(), text.size()));
}

TEST_CASE("CRC32 file helper matches memory") {
  const fs::path path = fs::temp_directory_path() / "covidx_crc32_test.bin";
  const std::string text = "123456789";
  {
    std::ofstream out(path, std::ios::binary);
    out.write(text.data(), static_cast<std::streamsize>(text.size()));
  }
  CHECK(covidx::detail::crc32_file(path.string()) == 0xCBF43926u);
  fs::remove(path);
}

TEST_CASE("CRC32 missing file throws IoError") {
  CHECK_THROWS_AS((void)covidx::detail::crc32_file("/nonexistent/covidx/file"),
                  covidx::IoError);
}

TEST_CASE("SHA-256 known vectors") {
  CHECK(covidx::detail::sha256_hex("abc", 3) ==
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
  CHECK(covidx::detail::sha256_hex("", 0) ==
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
}
