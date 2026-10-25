#include "doctest.h"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "covidx/detail/sha256.h"
#include "covidx/manifest.h"
#include "test_util.h"

namespace fs = std::filesystem;

namespace {

fs::path tiny_bam() {
  return fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam";
}
std::string build_bin() {
  return COVIDX_BUILD_BIN;
}
std::string query_bin() {
  return COVIDX_QUERY_BIN;
}
std::string quoted(const fs::path& path) {
  return "'" + path.string() + "'";
}

int run_build(const fs::path& out, const std::string& extra = "") {
  return tests::run_command_status(build_bin() + " --bam " +
                                   quoted(tiny_bam()) + " --out " +
                                   quoted(out) + " " + extra + " 2>/dev/null");
}

int run_query(const std::string& args, std::string* output) {
  return tests::run_command_status(query_bin() + " " + args, output);
}

struct Row {
  std::string contig;
  uint64_t start = 0;
  uint64_t end = 0;
  std::string unit;
  uint64_t features = 0;
  uint64_t covered = 0;
  uint64_t uncovered = 0;
  std::string fraction;
};

std::vector<Row> parse_rows(const std::string& text) {
  std::vector<Row> rows;
  std::istringstream in(text);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line.starts_with("#"))
      continue;
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (start <= line.size()) {
      const std::size_t tab = line.find('\t', start);
      fields.push_back(line.substr(
          start, tab == std::string::npos ? std::string::npos : tab - start));
      if (tab == std::string::npos)
        break;
      start = tab + 1;
    }
    REQUIRE(fields.size() == 8);
    rows.push_back(Row{fields[0], std::stoull(fields[1]),
                       std::stoull(fields[2]), fields[3],
                       std::stoull(fields[4]), std::stoull(fields[5]),
                       std::stoull(fields[6]), fields[7]});
  }
  return rows;
}

std::vector<bool> samtools_bitmap(const std::string& bam, uint64_t length) {
  std::string output;
  REQUIRE(tests::run_command_status(
              "samtools depth -aa -q 0 -Q 0 -G 0x800 -J -r chr1 " + bam,
              &output) == 0);
  std::vector<bool> covered(static_cast<std::size_t>(length), false);
  std::istringstream in(output);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (start <= line.size()) {
      const std::size_t tab = line.find('\t', start);
      fields.push_back(line.substr(
          start, tab == std::string::npos ? std::string::npos : tab - start));
      if (tab == std::string::npos)
        break;
      start = tab + 1;
    }
    REQUIRE(fields.size() >= 3);
    const uint64_t position = std::stoull(fields[1]) - 1;
    REQUIRE(position < length);
    covered[static_cast<std::size_t>(position)] = std::stoull(fields[2]) > 0;
  }
  return covered;
}

std::string blocks_bed(const std::vector<bool>& covered) {
  std::string text;
  std::size_t i = 0;
  while (i < covered.size()) {
    if (!covered[i]) {
      ++i;
      continue;
    }
    const std::size_t start = i;
    while (i < covered.size() && covered[i])
      ++i;
    text += "chr1\t" + std::to_string(start) + "\t" + std::to_string(i) + "\n";
  }
  return text;
}

} // namespace

TEST_CASE("CLI end-to-end bedtools coverage column 5 parity") {
  if (!tests::command_exists("samtools") ||
      !tests::command_exists("bedtools")) {
    MESSAGE("samtools/bedtools not found; skipping end-to-end parity");
    return;
  }
  tests::TempDir temp;
  const fs::path out = temp.path() / "index";
  REQUIRE(run_build(out, "--layout both") == 0);
  const fs::path index = out / "tiny.covidx";

  const uint64_t length = 1000;
  const std::vector<bool> covered =
      samtools_bitmap(tiny_bam().string(), length);
  const fs::path blocks = temp.path() / "blocks.bed";
  tests::write_file(blocks, blocks_bed(covered));

  std::vector<std::pair<uint64_t, uint64_t>> windows;
  for (uint64_t x = 0; x < length; ++x) {
    for (uint64_t len = 1; len <= 25 && x + len <= length; ++len) {
      windows.push_back({x, x + len});
    }
  }
  std::mt19937 rng(711);
  for (int i = 0; i < 500; ++i) {
    const uint64_t x = rng() % length;
    windows.push_back({x, x + 1 + rng() % (length - x)});
  }
  std::string windows_text;
  for (const auto& [x, y] : windows) {
    windows_text +=
        "chr1\t" + std::to_string(x) + "\t" + std::to_string(y) + "\n";
  }
  const fs::path windows_bed = temp.path() / "windows.bed";
  tests::write_file(windows_bed, windows_text);

  std::string bedtools_output;
  REQUIRE(tests::run_command_status("bedtools coverage -a " +
                                        quoted(windows_bed) + " -b " +
                                        quoted(blocks),
                                    &bedtools_output) == 0);
  std::vector<uint64_t> expected;
  std::istringstream bedtools_in(bedtools_output);
  std::string line;
  while (std::getline(bedtools_in, line)) {
    if (line.empty())
      continue;
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (start <= line.size()) {
      const std::size_t tab = line.find('\t', start);
      fields.push_back(line.substr(
          start, tab == std::string::npos ? std::string::npos : tab - start));
      if (tab == std::string::npos)
        break;
      start = tab + 1;
    }
    REQUIRE(fields.size() >= 7);
    expected.push_back(std::stoull(fields[4]));
  }
  REQUIRE(expected.size() == windows.size());

  std::string query_output;
  REQUIRE(run_query("--index " + quoted(index) + " --unit pooled --bed " +
                        quoted(windows_bed),
                    &query_output) == 0);
  const std::vector<Row> rows = parse_rows(query_output);
  REQUIRE(rows.size() == windows.size());
  for (std::size_t i = 0; i < rows.size(); ++i) {
    CHECK(rows[i].start == windows[i].first);
    CHECK(rows[i].end == windows[i].second);
    CHECK(rows[i].covered == expected[i]);
    CHECK(rows[i].uncovered == rows[i].end - rows[i].start - rows[i].covered);
  }
}

TEST_CASE("CLI scoped build matches a full build inside the scope") {
  tests::TempDir temp;
  const fs::path full_out = temp.path() / "full";
  const fs::path scoped_out = temp.path() / "scoped";
  REQUIRE(run_build(full_out) == 0);
  REQUIRE(run_build(scoped_out,
                    "--regions " + quoted(fs::path(COVIDX_TEST_DATA_DIR) /
                                          "tiny_scope.bed")) == 0);
  const fs::path bed = temp.path() / "windows.bed";
  tests::write_file(bed, "chr1\t100\t400\n"
                         "chr1\t400\t650\n"
                         "chr1\t900\t1000\n"
                         "chr2\t0\t500\n"
                         "chr1\t100\t300\n");
  std::string full_output;
  std::string scoped_output;
  REQUIRE(run_query("--index " + quoted(full_out / "tiny.covidx") +
                        " --unit S1 --bed " + quoted(bed),
                    &full_output) == 0);
  REQUIRE(run_query("--index " + quoted(scoped_out / "tiny.covidx") +
                        " --unit S1 --bed " + quoted(bed),
                    &scoped_output) == 0);
  const auto full_rows = parse_rows(full_output);
  const auto scoped_rows = parse_rows(scoped_output);
  REQUIRE(full_rows.size() == scoped_rows.size());
  for (std::size_t i = 0; i < full_rows.size(); ++i) {
    CHECK(full_rows[i].features == scoped_rows[i].features);
    CHECK(full_rows[i].covered == scoped_rows[i].covered);
  }
  const fs::path outside = temp.path() / "outside.bed";
  tests::write_file(outside, "chr1\t0\t1000\n");
  CHECK(run_query("--index " + quoted(scoped_out / "tiny.covidx") +
                      " --unit S1 --bed " + quoted(outside) + " 2>/dev/null",
                  nullptr) == 2);
}

TEST_CASE("README quickstart commands work") {
  tests::TempDir temp;
  const fs::path out = temp.path() / "indexes";
  CHECK(run_build(out) == 0);
  const fs::path index = out / "tiny.covidx";
  REQUIRE(fs::exists(index / "manifest.json"));

  std::string output;
  CHECK(
      run_query("--index " + quoted(index) + " --unit S1 --region chr1:1-1000",
                &output) == 0);
  CHECK(output.find("chr1\t0\t1000\tS1\t6\t60\t940\t0.060000") !=
        std::string::npos);

  output.clear();
  CHECK(run_query("--index " + quoted(index) + " --unit S1 --window 500",
                  &output) == 0);
  CHECK(output.find("chr1\t0\t500\tS1\t2\t20\t480\t0.040000") !=
        std::string::npos);

  const fs::path contigs_out = temp.path() / "contigs";
  CHECK(run_build(contigs_out, "--contigs chr1") == 0);
  CHECK(fs::exists(contigs_out / "tiny.covidx" / "manifest.json"));
}

TEST_CASE("format documentation matches the implemented constants") {
  const fs::path doc = fs::path(COVIDX_SOURCE_DIR) / "docs" / "format.md";
  REQUIRE(fs::exists(doc));
  const std::string text = tests::read_file(doc);
  CHECK(text.find("COVIDXU1") != std::string::npos);
  CHECK(text.find("48") != std::string::npos);
  CHECK(text.find("uint32") != std::string::npos);
  CHECK(text.find("R+1") != std::string::npos);
  CHECK(text.find("little-endian") != std::string::npos);
}
