#include "doctest.h"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "covidx/cov_index.h"
#include "covidx/error.h"
#include "covidx/manifest.h"
#include "test_util.h"

namespace fs = std::filesystem;
using covidx::CovIndex;
using covidx::OpenOptions;
using covidx::Validation;

namespace {

fs::path tiny_bam() {
  return fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam";
}

fs::path make_index(const fs::path& parent, const std::string& name = "tiny") {
  CovIndex index;
  index.build(tiny_bam().string());
  index.save(parent.string(), name);
  return parent / (name + ".covidx");
}

fs::path copy_index(const fs::path& source, const fs::path& parent) {
  fs::create_directories(parent);
  const fs::path target = parent / source.filename();
  fs::copy(source, target,
           fs::copy_options::recursive | fs::copy_options::overwrite_existing);
  return target;
}

std::string patch_byte(std::string bytes, std::size_t offset, uint8_t value) {
  REQUIRE(offset < bytes.size());
  bytes[offset] = static_cast<char>(value);
  return bytes;
}

std::string patch_u32(std::string bytes, std::size_t offset, uint32_t value) {
  REQUIRE(offset + 4 <= bytes.size());
  for (std::size_t i = 0; i < 4; ++i)
    bytes[offset + i] = static_cast<char>(value >> (8 * i));
  return bytes;
}

std::string patch_u64(std::string bytes, std::size_t offset, uint64_t value) {
  REQUIRE(offset + 8 <= bytes.size());
  for (std::size_t i = 0; i < 8; ++i)
    bytes[offset + i] = static_cast<char>(value >> (8 * i));
  return bytes;
}

void remove_sidecar(const fs::path& path) {
  fs::remove(path.string() + ".crc32");
}

void expect_format_error(const fs::path& index, const OpenOptions& options) {
  CAPTURE(index.string());
  CovIndex loaded;
  CHECK_THROWS_AS(loaded.load(index.string(), options), covidx::FormatError);
}

} // namespace

TEST_CASE("load rejects manifest and file-level corruption") {
  tests::TempDir temp;
  const fs::path index = make_index(temp.path());
  const fs::path unit = index / "data" / "u000000.bin";

  // Manifest CRC mismatch from a changed unit file.
  {
    const fs::path copy = copy_index(index, temp.path() / "crc");
    std::string bytes = tests::read_file(copy / "data" / "u000000.bin");
    bytes = patch_byte(bytes, 0, 'X');
    tests::write_file(copy / "data" / "u000000.bin", bytes);
    expect_format_error(copy, OpenOptions{});
  }

  // A disagreeing sidecar always fails.
  {
    const fs::path copy = copy_index(index, temp.path() / "sidecar");
    tests::write_file(copy / "data" / "u000000.bin.crc32", "00000000\n");
    expect_format_error(copy, OpenOptions{});
  }

  // A missing unit file fails.
  {
    const fs::path copy = copy_index(index, temp.path() / "missing");
    fs::remove(copy / "data" / "u000000.bin");
    expect_format_error(copy, OpenOptions{});
  }

  // A changed TSV fails after disabling CRC checks.
  {
    const fs::path copy = copy_index(index, temp.path() / "tsv");
    std::string text = tests::read_file(copy / "contigs.tsv");
    text[0] = text[0] == '0' ? '1' : '0';
    tests::write_file(copy / "contigs.tsv", text);
    remove_sidecar(copy / "contigs.tsv");
    OpenOptions options;
    options.validate_crc = false;
    expect_format_error(copy, options);
  }

  // Manifest scope corruption.
  {
    const fs::path copy = copy_index(index, temp.path() / "scope");
    std::string text = tests::read_file(copy / "manifest.json");
    const std::size_t pos = text.find("\"end\": 1000");
    REQUIRE(pos != std::string::npos);
    text.replace(pos, 11, "\"end\": 1001");
    tests::write_file(copy / "manifest.json", text);
    expect_format_error(copy, OpenOptions{});
  }

  // partial_scope corruption.
  {
    const fs::path copy = copy_index(index, temp.path() / "partial");
    std::string text = tests::read_file(copy / "manifest.json");
    const std::string old_value = "\"partial_scope\": false";
    const std::size_t pos = text.find(old_value);
    REQUIRE(pos != std::string::npos);
    text.replace(pos, old_value.size(), "\"partial_scope\": true");
    tests::write_file(copy / "manifest.json", text);
    expect_format_error(copy, OpenOptions{});
  }
}

TEST_CASE("unit-file corruption is caught by the right validation level") {
  tests::TempDir temp;
  const fs::path index = make_index(temp.path());
  const std::string good = tests::read_file(index / "data" / "u000000.bin");
  const std::size_t segment = 160; // 64-byte header + 2 * 48-byte metadata
  const std::size_t runs = 6;      // S1 on chr1

  OpenOptions no_crc;
  no_crc.validate_crc = false;
  no_crc.validation = Validation::Full;

  // Header corruption is caught even at Validation::None.
  {
    const fs::path copy = copy_index(index, temp.path() / "magic");
    tests::write_file(copy / "data" / "u000000.bin", patch_byte(good, 0, 'X'));
    remove_sidecar(copy / "data" / "u000000.bin");
    OpenOptions none;
    none.validate_crc = false;
    none.validation = Validation::None;
    expect_format_error(copy, none);
  }

  // Wrong contig count.
  {
    const fs::path copy = copy_index(index, temp.path() / "count");
    tests::write_file(copy / "data" / "u000000.bin", patch_u32(good, 12, 3));
    remove_sidecar(copy / "data" / "u000000.bin");
    expect_format_error(copy, no_crc);
  }

  // Unaligned segment offset.
  {
    const fs::path copy = copy_index(index, temp.path() / "align");
    tests::write_file(copy / "data" / "u000000.bin",
                      patch_u64(good, 64, segment + 1));
    remove_sidecar(copy / "data" / "u000000.bin");
    expect_format_error(copy, no_crc);
  }

  // A middle cum entry is not sampled by Structural but is caught by Full.
  {
    const fs::path copy = copy_index(index, temp.path() / "cum");
    const std::size_t cum_offset = segment + 2 * runs * 4;
    tests::write_file(copy / "data" / "u000000.bin",
                      patch_u64(good, cum_offset + 2 * 8, 999));
    remove_sidecar(copy / "data" / "u000000.bin");
    OpenOptions none;
    none.validate_crc = false;
    none.validation = Validation::None;
    CovIndex unvalidated;
    CHECK_NOTHROW(unvalidated.load(copy.string(), none));
    OpenOptions structural;
    structural.validate_crc = false;
    structural.validation = Validation::Structural;
    expect_format_error(copy, structural);
    OpenOptions full = structural;
    full.validation = Validation::Full;
    expect_format_error(copy, full);
  }

  // A middle bucket entry is only caught by Full.
  {
    const fs::path copy = copy_index(index, temp.path() / "bucket");
    const std::size_t first_offset = segment + 2 * runs * 4 + (runs + 1) * 8;
    tests::write_file(copy / "data" / "u000000.bin",
                      patch_u32(good, first_offset + 4, 0));
    remove_sidecar(copy / "data" / "u000000.bin");
    OpenOptions structural;
    structural.validate_crc = false;
    structural.validation = Validation::Structural;
    CovIndex validated;
    CHECK_NOTHROW(validated.load(copy.string(), structural));
    OpenOptions full = structural;
    full.validation = Validation::Full;
    expect_format_error(copy, full);
  }

  // A truncated unit file is caught by size validation.
  {
    const fs::path copy = copy_index(index, temp.path() / "truncated");
    fs::resize_file(copy / "data" / "u000000.bin",
                    fs::file_size(copy / "data" / "u000000.bin") - 8);
    expect_format_error(copy, OpenOptions{});
  }

  // A corrupted run endpoint is caught structurally.
  {
    const fs::path copy = copy_index(index, temp.path() / "endpoint");
    const std::size_t end_offset = segment + runs * 4;
    const uint64_t first_end = [&] {
      const auto bytes = tests::read_file(index / "data" / "u000000.bin");
      uint64_t value = 0;
      for (std::size_t i = 0; i < 8; ++i) {
        value |= static_cast<uint64_t>(
                     static_cast<unsigned char>(bytes[end_offset + i]))
                 << (8 * i);
      }
      return value;
    }();
    tests::write_file(
        copy / "data" / "u000000.bin",
        patch_u32(good, segment, static_cast<uint32_t>(first_end)));
    remove_sidecar(copy / "data" / "u000000.bin");
    expect_format_error(copy, OpenOptions{});
  }
}
