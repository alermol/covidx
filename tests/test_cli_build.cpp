#include "doctest.h"

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "covidx/cov_index.h"
#include "covidx/manifest.h"
#include "test_util.h"

namespace fs = std::filesystem;
using covidx::CovIndex;

namespace {

std::string build_bin() {
  return COVIDX_BUILD_BIN;
}
fs::path tiny_bam() {
  return fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam";
}

int run_build(const std::string& arguments, std::string* output = nullptr) {
  return tests::run_command_status(build_bin() + " " + arguments, output);
}

std::string quoted(const fs::path& path) {
  return "'" + path.string() + "'";
}

fs::path build_index(const fs::path& parent, const std::string& extra = "") {
  const int status = run_build("--bam " + quoted(tiny_bam()) + " --out " +
                               quoted(parent) + " " + extra);
  REQUIRE(status == 0);
  return parent / "tiny.covidx";
}

} // namespace

TEST_CASE("covidx-build help and version") {
  std::string output;
  CHECK(run_build("--help", &output) == 0);
  CHECK(output.find("--bam") != std::string::npos);
  CHECK(output.find("--regions") != std::string::npos);
  CHECK(output.find("--feature-counts") == std::string::npos);
  output.clear();
  CHECK(run_build("--version", &output) == 0);
  CHECK(output == "covidx-build 0.2.0 (format 1.1.0)\n");
}

TEST_CASE("covidx-build usage errors exit 1") {
  CHECK(run_build("") == 1);
  CHECK(run_build("--bam") == 1);
  CHECK(run_build("--bogus") == 1);
  CHECK(run_build("--bam /nope.bam --out /tmp/x --threads 0") == 1);
  CHECK(run_build("--bam /nope.bam --out /tmp/x --min-mapq 256") == 1);
  CHECK(run_build("--bam /nope.bam --out /tmp/x --exclude-flags 0x10000") == 1);
  CHECK(run_build("--bam /nope.bam --out /tmp/x --ambiguous-read-bases nope") ==
        1);
  CHECK(run_build("--bam /nope.bam --out /tmp/x --layout nope") == 1);
  CHECK(run_build("--bam /nope.bam --out /tmp/x --name bad/name") == 1);
  CHECK(run_build("--bam /nope.bam --out /tmp/x --group bad") == 1);
  CHECK(run_build("--bam= --out /tmp/x") == 1);
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out=") == 1);
}

TEST_CASE("covidx-build default build and metadata") {
  tests::TempDir temp;
  const fs::path directory = build_index(temp.path());
  REQUIRE(fs::exists(directory / "manifest.json"));
  CovIndex index;
  index.load(directory.string());
  CHECK(index.units().size() == 4);
  CHECK(index.covered_count("S1", "chr1", 0, 1000) == 60);
  CHECK(index.manifest().source.records_examined == tests::tiny_reads().size());
}

TEST_CASE("covidx-build always stores feature counts") {
  tests::TempDir temp;
  const fs::path out = temp.path() / "features";
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " + quoted(out) +
                  " --layout both --group trio=S1,S2") == 0);
  const fs::path feature_index = out / "tiny.covidx";
  const covidx::Manifest manifest =
      covidx::Manifest::load((feature_index / "manifest.json").string());
  CHECK(manifest.format_version == "1.1.0");
  REQUIRE(manifest.units.size() == 6);
  CHECK(manifest.units[0].total_features == 6);
  CHECK(manifest.units[4].name == "trio");
  CHECK(manifest.units[4].total_features == 10);
  CHECK(manifest.units[5].name == "pooled");
  CHECK(manifest.units[5].total_features == 12);

  CovIndex index;
  index.load(feature_index.string());
  CHECK(index.feature_count("S2", "chr1", 200, 212) == 1);

  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " +
                  quoted(temp.path() / "rejected") + " --feature-counts") == 1);
}

TEST_CASE("covidx-build layouts, groups, and groups-file") {
  tests::TempDir temp;
  const fs::path out = temp.path() / "out";
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " + quoted(out) +
                  " --layout both --group trio=S1,S2") == 0);
  CovIndex index;
  index.load((out / "tiny.covidx").string());
  REQUIRE(index.units().size() == 6);
  CHECK(index.unit("trio") != nullptr);
  CHECK(index.total_covered("trio", "chr1") == 103);

  const fs::path groups = temp.path() / "groups.tsv";
  tests::write_file(groups, "# comment\ntrio\tS1\tS2\n");
  const fs::path out2 = temp.path() / "out2";
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " + quoted(out2) +
                  " --layout both --groups-file " + quoted(groups)) == 0);
  CovIndex second;
  second.load((out2 / "tiny.covidx").string());
  CHECK(second.total_covered("trio", "chr1") == 103);
}

TEST_CASE("covidx-build scope filters") {
  tests::TempDir temp;
  const fs::path contigs_out = temp.path() / "contigs";
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " +
                  quoted(contigs_out) + " --contigs chr1") == 0);
  CovIndex contigs;
  contigs.load((contigs_out / "tiny.covidx").string());
  REQUIRE(contigs.contigs().size() == 1);
  CHECK(contigs.contig("chr1") != nullptr);
  CHECK(contigs.contig("chr2") == nullptr);
  CHECK_THROWS_AS((void)contigs.covered_count("S1", "chr2", 0, 1),
                  covidx::QueryError);

  const fs::path list = temp.path() / "contigs.txt";
  tests::write_file(list, "chr2\n");
  const fs::path union_out = temp.path() / "union";
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " +
                  quoted(union_out) + " --contigs chr1 --contigs-file " +
                  quoted(list)) == 0);
  CovIndex both;
  both.load((union_out / "tiny.covidx").string());
  CHECK(both.contigs().size() == 2);

  const fs::path empty = temp.path() / "empty.txt";
  tests::write_file(empty, "");
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " +
                  quoted(temp.path() / "x") + " --contigs-file " +
                  quoted(empty)) == 1);
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " +
                  quoted(temp.path() / "x") + " --contigs nope") == 1);

  const fs::path regions_out = temp.path() / "regions";
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " +
                  quoted(regions_out) + " --regions " +
                  quoted(fs::path(COVIDX_TEST_DATA_DIR) / "tiny_scope.bed")) ==
        0);
  CovIndex regions;
  regions.load((regions_out / "tiny.covidx").string());
  CHECK(regions.manifest().index.partial_scope);
  CHECK(regions.total_covered("S1", "chr1") == 40);
  CHECK_THROWS_AS((void)regions.covered_count("S1", "chr1", 0, 1000),
                  covidx::QueryError);

  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " +
                  quoted(temp.path() / "x") + " --contigs chr1 --regions " +
                  quoted(fs::path(COVIDX_TEST_DATA_DIR) / "tiny_scope.bed")) ==
        1);
}

TEST_CASE("covidx-build overwrite and output errors") {
  tests::TempDir temp;
  const fs::path out = temp.path() / "out";
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " + quoted(out)) ==
        0);
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " + quoted(out)) ==
        3);
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " + quoted(out) +
                  " --overwrite") == 0);

  const fs::path file = temp.path() / "not_a_directory";
  tests::write_file(file, "x");
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " + quoted(file)) ==
        3);
}

TEST_CASE("covidx-build input errors exit 2") {
  tests::TempDir temp;
  CHECK(run_build("--bam " + quoted(temp.path() / "missing.bam") + " --out " +
                  quoted(temp.path() / "out")) == 2);
  CHECK(run_build("--bam " +
                  quoted(fs::path(COVIDX_TEST_DATA_DIR) / "tiny_scope.bed") +
                  " --out " + quoted(temp.path() / "out")) == 2);
  const fs::path bad_bed = temp.path() / "bad.bed";
  tests::write_file(bad_bed, "chr1\t10\t10\n");
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " +
                  quoted(temp.path() / "out") + " --regions " +
                  quoted(bad_bed)) == 2);

  const fs::path unknown_bed = temp.path() / "unknown.bed";
  tests::write_file(unknown_bed, "nope\t0\t10\n");
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " +
                  quoted(temp.path() / "out") + " --regions " +
                  quoted(unknown_bed)) == 2);

  const fs::path beyond_bed = temp.path() / "beyond.bed";
  tests::write_file(beyond_bed, "chr1\t0\t1001\n");
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " +
                  quoted(temp.path() / "out") + " --regions " +
                  quoted(beyond_bed)) == 2);
}

TEST_CASE("covidx-build work directory, checksum, and verbose") {
  tests::TempDir temp;
  const fs::path out = temp.path() / "out";
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " + quoted(out)) ==
        0);
  bool work_left = false;
  if (fs::exists(out / ".covidx-work")) {
    for (const auto& entry : fs::directory_iterator(out / ".covidx-work")) {
      if (entry.is_directory() &&
          entry.path().filename().string().starts_with("covidx-work-")) {
        work_left = true;
      }
    }
  }
  CHECK_FALSE(work_left);

  const fs::path keep_out = temp.path() / "keep";
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " + quoted(keep_out) +
                  " --keep-work") == 0);
  bool kept = false;
  for (const auto& entry : fs::directory_iterator(keep_out / ".covidx-work")) {
    if (entry.is_directory() &&
        entry.path().filename().string().starts_with("covidx-work-")) {
      kept = true;
    }
  }
  CHECK(kept);

  const fs::path sums_out = temp.path() / "sums";
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " + quoted(sums_out) +
                  " --checksum-source") == 0);
  const covidx::Manifest manifest = covidx::Manifest::load(
      (sums_out / "tiny.covidx" / "manifest.json").string());
  REQUIRE(manifest.source.sha256.has_value());
  CHECK(manifest.source.sha256->size() == 64);

  std::string output;
  const fs::path verbose_out = temp.path() / "verbose";
  CHECK(run_build("--bam " + quoted(tiny_bam()) + " --out " +
                      quoted(verbose_out) + " --verbose 2>&1",
                  &output) == 0);
  CHECK(output.find("contigs=") != std::string::npos);
}

TEST_CASE("covidx-build refuses a stale index with --require-fresh-index") {
  tests::TempDir temp;
  const fs::path bam = temp.path() / "copy.bam";
  fs::copy_file(tiny_bam(), bam);
  fs::copy_file(tiny_bam().string() + ".bai", bam.string() + ".bai");
  const auto future = fs::file_time_type::clock::now() + std::chrono::hours(1);
  fs::last_write_time(bam, future);
  CHECK(run_build("--bam " + quoted(bam) + " --out " +
                  quoted(temp.path() / "out") + " --require-fresh-index") == 2);
  CHECK(run_build("--bam " + quoted(bam) + " --out " +
                  quoted(temp.path() / "out2")) == 0);
}
