#include "doctest.h"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <random>
#include <string>
#include <vector>

#include "covidx/detail/bam_reader.h"
#include "covidx/detail/build.h"
#include "covidx/detail/run_index.h"
#include "test_util.h"

namespace fs = std::filesystem;
using covidx::BuildOptions;
using covidx::detail::BamHeaderInfo;
using covidx::detail::RunIndex;

namespace {

std::vector<bool> samtools_bitmap(const std::string& bam,
                                  const std::string& contig, uint64_t length,
                                  const std::string& flags) {
  std::string output;
  const std::string command =
      "samtools depth -aa -q 0 " + flags + " -r " + contig + " " + bam;
  REQUIRE(tests::run_command(command, &output) == 0);
  std::vector<bool> covered(static_cast<std::size_t>(length), false);
  std::size_t start = 0;
  while (start < output.size()) {
    const std::size_t end = output.find('\n', start);
    const std::string line = output.substr(
        start, end == std::string::npos ? std::string::npos : end - start);
    start = end == std::string::npos ? output.size() : end + 1;
    if (line.empty())
      continue;
    const std::vector<std::string> fields = tests::split_tabs(line);
    const uint64_t pos = std::stoull(fields[1]) - 1;
    const uint64_t depth = std::stoull(fields[2]);
    REQUIRE(pos < length);
    covered[static_cast<std::size_t>(pos)] = depth > 0;
  }
  return covered;
}

std::string blocks_bed(const std::vector<bool>& covered) {
  std::string out;
  std::size_t i = 0;
  while (i < covered.size()) {
    if (!covered[i]) {
      ++i;
      continue;
    }
    const std::size_t start = i;
    while (i < covered.size() && covered[i])
      ++i;
    out += "chr1\t" + std::to_string(start) + "\t" + std::to_string(i) + "\n";
  }
  return out;
}

struct Window {
  uint64_t x;
  uint64_t y;
};

std::vector<Window> make_windows(uint64_t length) {
  std::vector<Window> windows;
  for (uint64_t x = 0; x < length; ++x) {
    for (uint64_t len = 1; len <= 20 && x + len <= length; ++len) {
      windows.push_back(Window{x, x + len});
    }
  }
  std::mt19937 rng(9001);
  for (int i = 0; i < 1000; ++i) {
    const uint64_t x = rng() % length;
    const uint64_t y = x + 1 + rng() % (length - x);
    windows.push_back(Window{x, y});
  }
  return windows;
}

void check_bedtools(const fs::path& bam_path, const std::string& flags,
                    const BuildOptions& options) {
  if (!tests::command_exists("bedtools")) {
    MESSAGE("bedtools not found; skipping parity test");
    return;
  }
  const BamHeaderInfo header =
      covidx::detail::read_bam_header(bam_path.string());
  const uint64_t length = header.contigs[0].length;
  const std::vector<bool> covered =
      samtools_bitmap(bam_path.string(), "chr1", length, flags);
  RunIndex index;
  covidx::detail::build_single_unit(header, options, 0, index);

  tests::TempDir temp;
  const fs::path blocks = temp.path() / "blocks.bed";
  tests::write_file(blocks, blocks_bed(covered));

  const std::vector<Window> windows = make_windows(length);
  std::string windows_text;
  for (const Window& window : windows) {
    windows_text += "chr1\t" + std::to_string(window.x) + "\t" +
                    std::to_string(window.y) + "\n";
  }
  const fs::path windows_bed = temp.path() / "windows.bed";
  tests::write_file(windows_bed, windows_text);

  const std::string command =
      "bedtools coverage -a " + windows_bed.string() + " -b " + blocks.string();
  std::string output;
  REQUIRE(tests::run_command(command, &output) == 0);

  std::vector<std::vector<std::string>> rows;
  std::size_t start = 0;
  while (start < output.size()) {
    const std::size_t end = output.find('\n', start);
    const std::string line = output.substr(
        start, end == std::string::npos ? std::string::npos : end - start);
    start = end == std::string::npos ? output.size() : end + 1;
    if (!line.empty())
      rows.push_back(tests::split_tabs(line));
  }
  REQUIRE(rows.size() == windows.size());
  for (std::size_t i = 0; i < windows.size(); ++i) {
    const auto& row = rows[i];
    REQUIRE(row.size() >= 7);
    const uint64_t got = std::stoull(row[4]); // num_bases_covered
    const uint64_t length_field = std::stoull(row[5]);
    const uint64_t want = index.covered_count(windows[i].x, windows[i].y);
    CHECK(length_field == windows[i].y - windows[i].x);
    CHECK(got == want);
  }
}

} // namespace

TEST_CASE("bedtools coverage column 5 parity (default policy)") {
  if (!tests::command_exists("samtools")) {
    MESSAGE("samtools not found; skipping parity test");
    return;
  }
  const fs::path bam = fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam";
  check_bedtools(bam, "-Q 0 -G 0x800 -J", BuildOptions{});
}

TEST_CASE("bedtools coverage column 5 parity (deletions are gaps)") {
  if (!tests::command_exists("samtools")) {
    MESSAGE("samtools not found; skipping parity test");
    return;
  }
  const fs::path bam = fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam";
  BuildOptions options;
  options.count_deletions = false;
  check_bedtools(bam, "-Q 0 -G 0x800", options);
}
