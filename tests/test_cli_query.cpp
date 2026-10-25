#include "doctest.h"

#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "covidx/cov_index.h"
#include "test_util.h"

namespace fs = std::filesystem;
using covidx::CovIndex;

namespace {

std::string query_bin() {
  return COVIDX_QUERY_BIN;
}
std::string build_bin() {
  return COVIDX_BUILD_BIN;
}
fs::path tiny_bam() {
  return fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam";
}

std::string quoted(const fs::path& path) {
  return "'" + path.string() + "'";
}

int run_query(const std::string& arguments, std::string* output = nullptr) {
  return tests::run_command_status(query_bin() + " " + arguments, output);
}

int run_build(const std::string& arguments, std::string* output = nullptr) {
  return tests::run_command_status(build_bin() + " " + arguments, output);
}

fs::path make_tiny_index(const fs::path& parent,
                         const std::string& extra = "") {
  const int status = run_build("--bam " + quoted(tiny_bam()) + " --out " +
                               quoted(parent) + " " + extra);
  if (status != 0)
    throw std::runtime_error("covidx-build failed with status " +
                             std::to_string(status));
  return parent / "tiny.covidx";
}

fs::path make_custom_bam(const fs::path& parent, std::string_view header,
                         std::vector<tests::ReadSpec> reads,
                         const std::string& name = "custom") {
  const fs::path bam = parent / (name + ".bam");
  tests::make_bam(bam, header, std::move(reads));
  const fs::path out = parent / (name + "-index");
  const int status =
      run_build("--bam " + quoted(bam) + " --out " + quoted(out));
  if (status != 0)
    throw std::runtime_error("custom build failed");
  return out / (name + ".covidx");
}

constexpr std::string_view kColonHeader = "@HD\tVN:1.6\tSO:coordinate\n"
                                          "@SQ\tSN:chr:1\tLN:100\n"
                                          "@RG\tID:rg1\tSM:S1\n";

constexpr std::string_view kAmbiguousHeader = "@HD\tVN:1.6\tSO:coordinate\n"
                                              "@SQ\tSN:chr1\tLN:100\n"
                                              "@SQ\tSN:chr1:1-10\tLN:50\n"
                                              "@RG\tID:rg1\tSM:S1\n";

} // namespace

TEST_CASE("covidx-query help and version") {
  std::string output;
  CHECK(run_query("--help", &output) == 0);
  CHECK(output.find("--region") != std::string::npos);
  CHECK(output.find("--window") != std::string::npos);
  CHECK(output.find("--features") == std::string::npos);
  output.clear();
  CHECK(run_query("--version", &output) == 0);
  CHECK(output == "covidx-query 0.2.0 (format 1.1.0)\n");
}

TEST_CASE("covidx-query list modes") {
  tests::TempDir temp;
  const fs::path index = make_tiny_index(temp.path());
  std::string output;
  CHECK(run_query("--index " + quoted(index) + " --list-units", &output) == 0);
  CHECK(output ==
        "# unit\tkind\nS1\tsample\nS2\tsample\nS3\tsample\nS4\tsample\n");
  output.clear();
  CHECK(run_query("--index " + quoted(index) + " --list-contigs", &output) ==
        0);
  CHECK(output == "# contig\tlength\nchr1\t1000\nchr2\t500\n");
  output.clear();
  CHECK(run_query("--index " + quoted(index) + " --list-units --no-header",
                  &output) == 0);
  CHECK(output == "S1\tsample\nS2\tsample\nS3\tsample\nS4\tsample\n");
}

TEST_CASE("covidx-query exact region and tiling output") {
  tests::TempDir temp;
  const fs::path index = make_tiny_index(temp.path());
  std::string output;
  CHECK(run_query("--index " + quoted(index) + " --unit S1 --region chr1:1-10",
                  &output) == 0);
  CHECK(output ==
        "#chrom\tstart\tend\tunit\tnum_features\tcovered\tuncovered\t"
        "covered_fraction\n"
        "chr1\t0\t10\tS1\t1\t10\t0\t1.000000\n");

  output.clear();
  CHECK(run_query("--index " + quoted(index) +
                      " --unit S1 --region chr1:1-1000 --window 400 --step 400",
                  &output) == 0);
  CHECK(output ==
        "#chrom\tstart\tend\tunit\tnum_features\tcovered\tuncovered\t"
        "covered_fraction\n"
        "chr1\t0\t400\tS1\t1\t10\t390\t0.025000\n"
        "chr1\t400\t800\tS1\t3\t30\t370\t0.075000\n"
        "chr1\t800\t1000\tS1\t2\t20\t180\t0.100000\n");

  output.clear();
  CHECK(run_query("--index " + quoted(index) + " --unit S1 --window 400",
                  &output) == 0);
  CHECK(output ==
        "#chrom\tstart\tend\tunit\tnum_features\tcovered\tuncovered\t"
        "covered_fraction\n"
        "chr1\t0\t400\tS1\t1\t10\t390\t0.025000\n"
        "chr1\t400\t800\tS1\t3\t30\t370\t0.075000\n"
        "chr1\t800\t1000\tS1\t2\t20\t180\t0.100000\n"
        "chr2\t0\t400\tS1\t0\t0\t400\t0.000000\n"
        "chr2\t400\t500\tS1\t0\t0\t100\t0.000000\n");

  output.clear();
  CHECK(run_query("--index " + quoted(index) +
                      " --unit S1 --window 400 --step 600 --contig chr1",
                  &output) == 0);
  CHECK(output ==
        "#chrom\tstart\tend\tunit\tnum_features\tcovered\tuncovered\t"
        "covered_fraction\n"
        "chr1\t0\t400\tS1\t1\t10\t390\t0.025000\n"
        "chr1\t600\t1000\tS1\t3\t30\t370\t0.075000\n");

  // A window wider than the region must clamp, not overflow the tiling
  // arithmetic. The region starts at a non-zero coordinate so the old
  // position + window expression wrapped.
  output.clear();
  CHECK(
      run_query(
          "--index " + quoted(index) +
              " --unit S1 --region chr1:101-1000 --window 18446744073709551615",
          &output) == 0);
  CHECK(output ==
        "#chrom\tstart\tend\tunit\tnum_features\tcovered\tuncovered\t"
        "covered_fraction\n"
        "chr1\t100\t1000\tS1\t5\t50\t850\t0.055556\n");
}

TEST_CASE("covidx-query reads BED before truncating the output path") {
  tests::TempDir temp;
  const fs::path index = make_tiny_index(temp.path());
  const fs::path bed = temp.path() / "same-path.bed";
  tests::write_file(bed, "chr1\t0\t10\n");
  CHECK(run_query("--index " + quoted(index) + " --unit S1 --bed " +
                  quoted(bed) + " --output " + quoted(bed)) == 0);
  const std::string output = tests::read_file(bed);
  CHECK(output.find("chr1\t0\t10\tS1\t1\t10\t0\t1.000000\n") !=
        std::string::npos);
}

TEST_CASE("covidx-query region conventions and errors") {
  tests::TempDir temp;
  const fs::path index = make_tiny_index(temp.path());
  std::string output;
  CHECK(run_query("--index " + quoted(index) +
                      " --unit S1 --zero-based --region chr1:0-10",
                  &output) == 0);
  CHECK(output.find("chr1\t0\t10\tS1\t1\t10\t0\t1.000000") !=
        std::string::npos);
  output.clear();
  CHECK(run_query("--index " + quoted(index) +
                      " --unit S1 --one-based --region chr1:1-10",
                  &output) == 0);
  CHECK(output.find("chr1\t0\t10\tS1\t1\t10\t0\t1.000000") !=
        std::string::npos);
  CHECK(run_query("--index " + quoted(index) +
                  " --unit S1 --region chr1:0-10 2>/dev/null") == 1);
  CHECK(run_query("--index " + quoted(index) +
                  " --unit S1 --region chr1:11-10 2>/dev/null") == 1);
  CHECK(run_query("--index " + quoted(index) +
                  " --unit S1 --region chr1:1-1001 2>/dev/null") == 1);
  CHECK(run_query("--index " + quoted(index) +
                  " --unit S1 --region nope 2>/dev/null") == 2);
  CHECK(run_query("--index " + quoted(index) +
                  " --unit S1 --region chr1:abc 2>/dev/null") == 2);
  CHECK(run_query("--index " + quoted(index) +
                  " --unit S1 --region chr1 2>/dev/null") == 0);
  CHECK(run_query("--index " + quoted(index) +
                  " --unit S1 --region chr1 --zero-based "
                  "--one-based 2>/dev/null") == 1);
}

TEST_CASE("covidx-query BED parsing and validation") {
  tests::TempDir temp;
  const fs::path index = make_tiny_index(temp.path());
  const fs::path bed = temp.path() / "query.bed";
  tests::write_file(bed, "# comment\n"
                         "track name=x\n"
                         "browser position chr1\n"
                         "\n"
                         "chr1\t0\t10\tname\textra\r\n"
                         "chr1\t10\t20\n"
                         "chr1\t0\t1000\n");
  std::string output;
  CHECK(
      run_query("--index " + quoted(index) + " --unit S1 --bed " + quoted(bed),
                &output) == 0);
  CHECK(output ==
        "#chrom\tstart\tend\tunit\tnum_features\tcovered\tuncovered\t"
        "covered_fraction\n"
        "chr1\t0\t10\tS1\t1\t10\t0\t1.000000\n"
        "chr1\t10\t20\tS1\t0\t0\t10\t0.000000\n"
        "chr1\t0\t1000\tS1\t6\t60\t940\t0.060000\n");

  output.clear();
  const std::string stdin_command = "printf 'chr1\\t0\\t10\\n' | " +
                                    query_bin() + " --index " + quoted(index) +
                                    " --unit S1 --bed -";
  CHECK(tests::run_command_status(stdin_command, &output) == 0);
  CHECK(output.find("chr1\t0\t10\tS1\t1\t10\t0\t1.000000") !=
        std::string::npos);

  const fs::path unknown = temp.path() / "unknown.bed";
  tests::write_file(unknown, "nope\t0\t10\n");
  CHECK(run_query("--index " + quoted(index) + " --unit S1 --bed " +
                  quoted(unknown) + " 2>/dev/null") == 2);
  const fs::path out_of_range = temp.path() / "range.bed";
  tests::write_file(out_of_range, "chr1\t0\t1001\n");
  CHECK(run_query("--index " + quoted(index) + " --unit S1 --bed " +
                  quoted(out_of_range) + " 2>/dev/null") == 2);
  const fs::path zero = temp.path() / "zero.bed";
  tests::write_file(zero, "chr1\t10\t10\n");
  CHECK(run_query("--index " + quoted(index) + " --unit S1 --bed " +
                  quoted(zero) + " 2>/dev/null") == 2);
  const fs::path bad = temp.path() / "bad.bed";
  tests::write_file(bad, "chr1\tx\t10\n");
  CHECK(run_query("--index " + quoted(index) + " --unit S1 --bed " +
                  quoted(bad) + " 2>/dev/null") == 2);
}

TEST_CASE("covidx-query unit resolution and output options") {
  tests::TempDir temp;
  const fs::path index = make_tiny_index(temp.path());
  CHECK(run_query("--index " + quoted(index) + " --window 100 2>/dev/null") ==
        1);
  std::string output;
  CHECK(run_query("--index " + quoted(index) +
                      " --unit S1 --window 1000 --no-header",
                  &output) == 0);
  CHECK(output == "chr1\t0\t1000\tS1\t6\t60\t940\t0."
                  "060000\nchr2\t0\t500\tS1\t0\t0\t500\t0.000000\n");

  const fs::path out = temp.path() / "result.tsv";
  CHECK(run_query("--index " + quoted(index) +
                  " --unit S1 --window 1000 --output " + quoted(out)) == 0);
  CHECK(tests::read_file(out) ==
        "#chrom\tstart\tend\tunit\tnum_features\tcovered\tuncovered\t"
        "covered_fraction\n"
        "chr1\t0\t1000\tS1\t6\t60\t940\t0.060000\n"
        "chr2\t0\t500\tS1\t0\t0\t500\t0.000000\n");

  const fs::path pooled_out = temp.path() / "pooled";
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " +
                  quoted(pooled_out) + " --layout pooled") == 0);
  output.clear();
  CHECK(run_query("--index " + quoted(pooled_out / "tiny.covidx") +
                      " --window 1000",
                  &output) == 0);
  CHECK(output.find("\tpooled\t12\t123\t") != std::string::npos);
}

TEST_CASE("covidx-query zero-unit index exits 2") {
  tests::TempDir temp;
  const std::string header =
      "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n";
  const fs::path bam = temp.path() / "empty.bam";
  tests::make_bam(bam, header, {});
  const fs::path out = temp.path() / "out";
  CHECK(run_build("--bam " + quoted(bam) + " --out " + quoted(out)) == 0);
  CHECK(run_query("--index " + quoted(out / "empty.covidx") +
                  " --window 10 2>/dev/null") == 2);
}

TEST_CASE("covidx-query partial scope tiling and errors") {
  tests::TempDir temp;
  const fs::path index = make_tiny_index(
      temp.path(),
      "--regions " + quoted(fs::path(COVIDX_TEST_DATA_DIR) / "tiny_scope.bed"));
  std::string output;
  CHECK(run_query("--index " + quoted(index) + " --unit S1 --window 300",
                  &output) == 0);
  CHECK(output ==
        "#chrom\tstart\tend\tunit\tnum_features\tcovered\tuncovered\t"
        "covered_fraction\n"
        "chr1\t100\t400\tS1\t0\t0\t300\t0.000000\n"
        "chr1\t400\t650\tS1\t3\t30\t220\t0.120000\n"
        "chr1\t900\t1000\tS1\t1\t10\t90\t0.100000\n"
        "chr2\t0\t300\tS1\t0\t0\t300\t0.000000\n"
        "chr2\t300\t500\tS1\t0\t0\t200\t0.000000\n");
  output.clear();
  CHECK(run_query("--index " + quoted(index) +
                      " --unit S1 --region chr1:101-650 --window 300",
                  &output) == 0);
  CHECK(output ==
        "#chrom\tstart\tend\tunit\tnum_features\tcovered\tuncovered\t"
        "covered_fraction\n"
        "chr1\t100\t400\tS1\t0\t0\t300\t0.000000\n"
        "chr1\t400\t650\tS1\t3\t30\t220\t0.120000\n");
  CHECK(run_query("--index " + quoted(index) +
                  " --unit S1 --region chr1 2>/dev/null") == 2);
  CHECK(run_query("--index " + quoted(index) +
                  " --unit S1 --region chr1:700-800 2>/dev/null") == 2);
  const fs::path bed = temp.path() / "scope.bed";
  tests::write_file(bed, "chr1\t0\t1000\n");
  CHECK(run_query("--index " + quoted(index) + " --unit S1 --bed " +
                  quoted(bed) + " 2>/dev/null") == 2);
}

TEST_CASE("covidx-query contig names containing colons") {
  tests::TempDir temp;
  tests::ReadSpec read{
      "r1", 0, 0, 0, 60, {{'M', 10}}, tests::repeat_sequence(10), "rg1"};
  const fs::path index =
      make_custom_bam(temp.path(), kColonHeader, {read}, "colon");
  std::string output;
  CHECK(run_query("--index " + quoted(index) + " --unit S1 --region 'chr:1'",
                  &output) == 0);
  CHECK(output.find("chr:1\t0\t100\tS1\t1\t10\t90\t0.100000") !=
        std::string::npos);
  output.clear();
  CHECK(
      run_query("--index " + quoted(index) + " --unit S1 --region 'chr:1:1-10'",
                &output) == 0);
  CHECK(output.find("chr:1\t0\t10\tS1\t1\t10\t0\t1.000000") !=
        std::string::npos);

  tests::ReadSpec read1 = read;
  read1.tid = 0;
  tests::ReadSpec read2 = read;
  read2.tid = 1;
  const fs::path ambiguous = make_custom_bam(temp.path(), kAmbiguousHeader,
                                             {read1, read2}, "ambiguous");
  CHECK(run_query("--index " + quoted(ambiguous) +
                  " --unit S1 --region 'chr1:1-10' "
                  "2>/dev/null") == 2);
}

TEST_CASE("covidx-query option validation exits 1") {
  tests::TempDir temp;
  const fs::path index = make_tiny_index(temp.path());
  CHECK(run_query("--index= --unit S1 --region chr1:1-10") == 1);
  const std::string base = "--index " + quoted(index);
  CHECK(run_query("--index " + quoted(index) +
                  " --one-based --zero-based --region chr1:1-10 "
                  "2>/dev/null") == 1);
  CHECK(run_query(base + " --bed x --window 10 2>/dev/null") == 1);
  CHECK(run_query(base + " --step 10 --window 0 2>/dev/null") == 1);
  CHECK(run_query(
            base +
            " --window 10 --contig chr1 --region chr1:1-10 2>/dev/null") == 1);
  CHECK(run_query(base + " --list-units --window 10 2>/dev/null") == 1);
  CHECK(run_query(base + " --window 10 --features 2>/dev/null") == 1);
  CHECK(run_query(base + " --validation nope --window 10 2>/dev/null") == 1);
  CHECK(run_query(base + " --unit S1 --window 10 --contig nope 2>/dev/null") ==
        2);
  CHECK(run_query("--index " + quoted(temp.path() / "missing") +
                  " --window 10 2>/dev/null") == 2);
}

TEST_CASE("covidx-query always outputs feature counts") {
  tests::TempDir temp;
  const fs::path features =
      make_tiny_index(temp.path() / "features",
                      "--layout both --group trio=S1,S2");
  std::string output;
  CHECK(run_query("--index " + quoted(features) +
                      " --unit S1 --region chr1:1-1000",
                  &output) == 0);
  CHECK(output == "#chrom\tstart\tend\tunit\tnum_features\tcovered\tuncovered\t"
                  "covered_fraction\n"
                  "chr1\t0\t1000\tS1\t6\t60\t940\t0.060000\n");

  output.clear();
  CHECK(run_query("--index " + quoted(features) +
                      " --unit trio --region chr1:1-1000",
                  &output) == 0);
  CHECK(output.find("chr1\t0\t1000\ttrio\t10\t103\t897\t0.103000\n") !=
        std::string::npos);

  output.clear();
  CHECK(
      run_query("--index " + quoted(features) + " --unit S1 --region chr1:1-10",
                &output) == 0);
  CHECK(output ==
        "#chrom\tstart\tend\tunit\tnum_features\tcovered\tuncovered\t"
        "covered_fraction\n"
        "chr1\t0\t10\tS1\t1\t10\t0\t1.000000\n");
}

TEST_CASE("covidx-query loading options") {
  tests::TempDir temp;
  const fs::path index = make_tiny_index(temp.path());
  std::string output;
  CHECK(run_query("--index " + quoted(index) +
                      " --unit S1 --window 1000 --no-mmap",
                  &output) == 0);
  CHECK(output.find("chr1\t0\t1000\tS1\t6\t60\t") != std::string::npos);
  output.clear();
  CHECK(run_query("--index " + quoted(index) +
                      " --unit S1 --window 1000 --no-validate-crc "
                      "--validation full",
                  &output) == 0);
  CHECK(output.find("chr1\t0\t1000\tS1\t6\t60\t") != std::string::npos);
}
