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

struct PolicyCase {
  const char* name;
  BuildOptions options;
  std::string samtools_flags;
};

std::vector<bool> samtools_bitmap(const std::string& bam,
                                  const std::string& contig, uint64_t length,
                                  const std::string& flags) {
  const std::string command =
      "samtools depth -aa -q 0 " + flags + " -r " + contig + " " + bam;
  std::string output;
  const int status = tests::run_command(command, &output);
  REQUIRE(status == 0);
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
    REQUIRE(fields.size() >= 3);
    const uint64_t pos = std::stoull(fields[1]) - 1;
    const uint64_t depth = std::stoull(fields[2]);
    REQUIRE(pos < length);
    covered[static_cast<std::size_t>(pos)] = depth > 0;
  }
  return covered;
}

void check_against_bitmap(const RunIndex& index,
                          const std::vector<bool>& covered) {
  uint64_t total = 0;
  for (bool value : covered)
    total += value ? 1 : 0;
  CHECK(index.total_covered() == total);
  for (uint64_t x = 0; x < covered.size(); ++x) {
    uint64_t running = 0;
    const uint64_t limit = std::min<uint64_t>(covered.size(), x + 80);
    for (uint64_t y = x + 1; y <= limit; ++y) {
      running += covered[static_cast<std::size_t>(y - 1)] ? 1 : 0;
      CHECK(index.covered_count(x, y) == running);
    }
  }
  std::mt19937 rng(5150);
  for (int i = 0; i < 500; ++i) {
    const uint64_t x = rng() % covered.size();
    const uint64_t y = x + 1 + rng() % (covered.size() - x);
    uint64_t want = 0;
    for (uint64_t p = x; p < y; ++p)
      want += covered[static_cast<std::size_t>(p)] ? 1 : 0;
    CHECK(index.covered_count(x, y) == want);
  }
}

} // namespace

TEST_CASE("samtools depth parity across policy rows") {
  if (!tests::command_exists("samtools")) {
    MESSAGE("samtools not found; skipping parity test");
    return;
  }
  const fs::path bam = fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam";
  const BamHeaderInfo header = covidx::detail::read_bam_header(bam.string());

  std::vector<PolicyCase> cases;
  cases.push_back({"defaults", BuildOptions{}, "-Q 0 -G 0x800 -J"});
  {
    BuildOptions options;
    options.include_supplementary = true;
    cases.push_back({"include supplementary", options, "-Q 0 -J"});
  }
  {
    BuildOptions options;
    options.count_deletions = false;
    cases.push_back({"no deletions", options, "-Q 0 -G 0x800"});
  }
  {
    BuildOptions options;
    options.min_mapq = 20;
    cases.push_back({"min mapq 20", options, "-Q 20 -G 0x800 -J"});
  }
  {
    BuildOptions options;
    options.exclude_flags = 0;
    cases.push_back(
        {"no flag filtering", options, "-Q 0 -g 0x704 -G 0x800 -J"});
  }
  {
    BuildOptions options;
    options.exclude_flags = 0x604;
    cases.push_back(
        {"include secondary", options, "-Q 0 -g 0x100 -G 0x800 -J"});
  }

  for (const PolicyCase& test : cases) {
    CAPTURE(test.name);
    const std::vector<bool> expected = samtools_bitmap(
        bam.string(), "chr1", header.contigs[0].length, test.samtools_flags);
    RunIndex index;
    covidx::detail::build_single_unit(header, test.options, 0, index);
    check_against_bitmap(index, expected);
  }
}
