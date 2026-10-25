#include "doctest.h"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "covidx/cov_index.h"
#include "covidx/error.h"
#include "covidx/manifest.h"
#include "test_util.h"

namespace fs = std::filesystem;
using covidx::BuildOptions;
using covidx::CovIndex;
using covidx::Layout;
using covidx::OpenOptions;
using covidx::Region;

namespace {

fs::path tiny_bam() {
  return fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam";
}

fs::path find_work_dir(const fs::path& parent) {
  for (const auto& entry : fs::directory_iterator(parent)) {
    if (entry.is_directory() &&
        entry.path().filename().string().starts_with("covidx-work-")) {
      return entry.path();
    }
  }
  return {};
}

CovIndex build_tiny(const BuildOptions& options = {}) {
  CovIndex index;
  index.build(tiny_bam().string(), options);
  return index;
}

} // namespace

TEST_CASE("CovIndex empty state") {
  CovIndex index;
  CHECK_THROWS_AS((void)index.covered_count("S1", "chr1", 0, 10),
                  covidx::QueryError);
  CHECK_THROWS_AS((void)index.covered_fraction("S1", "chr1", 0, 10),
                  covidx::QueryError);
  CHECK_THROWS_AS((void)index.total_covered("S1", "chr1"), covidx::QueryError);
  CHECK_THROWS_AS((void)index.units(), covidx::QueryError);
  CHECK_THROWS_AS((void)index.contigs(), covidx::QueryError);
  CHECK_THROWS_AS((void)index.scope("chr1"), covidx::QueryError);
  CHECK_THROWS_AS((void)index.manifest(), covidx::QueryError);
  CHECK(index.unit("S1") == nullptr);
  CHECK(index.contig("chr1") == nullptr);
  CHECK(index.build_stats().records_examined == 0);
}

TEST_CASE("CovIndex build exposes metadata and exact counts") {
  const CovIndex index = build_tiny();
  REQUIRE(index.units().size() == 4);
  CHECK(index.units()[0].name == "S1");
  CHECK(index.units()[2].name == "S3");
  REQUIRE(index.contigs().size() == 2);
  CHECK(index.contigs()[0].name == "chr1");
  CHECK(index.contigs()[0].length == 1000);
  CHECK(index.unit("S3") != nullptr);
  CHECK(index.unit("S3")->total_covered == 0);
  CHECK(index.contig("chr2") != nullptr);
  CHECK(index.contig("nope") == nullptr);

  CHECK(index.total_covered("S1", "chr1") == 60);
  CHECK(index.total_covered("S2", "chr1") == 43);
  CHECK(index.total_covered("S3", "chr1") == 0);
  CHECK(index.total_covered("S4", "chr1") == 20);
  CHECK(index.total_covered("S1", "chr2") == 0);
  CHECK(index.covered_count("S1", "chr1", 0, 10) == 10);
  CHECK(index.covered_count("S1", "chr1", 10, 400) == 0);
  CHECK(index.covered_count("S1", "chr1", 0, 1000) == 60);
  CHECK(index.covered_count("S2", "chr1", 200, 212) == 12);
  CHECK(index.covered_count("S4", "chr1", 960, 970) == 10);
  CHECK(index.covered_fraction("S1", "chr1", 0, 1000) == doctest::Approx(0.06));
  CHECK(index.covered_fraction("S3", "chr1", 0, 1000) == 0.0);
  CHECK(index.covered_fraction("S1", "chr1", 0, 10) == 1.0);

  CHECK(index.manifest().format == "covidx");
  CHECK(index.manifest().layout == "per-sample");
  CHECK(index.manifest().policy.ambiguous_read_bases == "cover");
  CHECK(index.manifest().policy.effective_exclude_flags == 0x0F04);
  CHECK(index.manifest().index.total_runs == 13);
  CHECK(index.build_stats().records_examined == tests::tiny_reads().size());
  CHECK(index.build_stats().records_accepted == 12);
}

TEST_CASE("CovIndex query validation errors") {
  const CovIndex index = build_tiny();
  CHECK_THROWS_AS((void)index.covered_count("nope", "chr1", 0, 10),
                  covidx::QueryError);
  CHECK_THROWS_AS((void)index.covered_count("S1", "nope", 0, 10),
                  covidx::QueryError);
  CHECK_THROWS_AS((void)index.covered_count("S1", "chr1", 5, 5),
                  covidx::QueryError);
  CHECK_THROWS_AS((void)index.covered_count("S1", "chr1", 6, 5),
                  covidx::QueryError);
  CHECK_THROWS_AS((void)index.covered_count("S1", "chr1", 0, 1001),
                  covidx::QueryError);
  CHECK_THROWS_AS((void)index.total_covered("nope", "chr1"),
                  covidx::QueryError);
  CHECK_THROWS_AS((void)index.total_covered("S1", "nope"), covidx::QueryError);
  CHECK_THROWS_AS((void)index.scope("nope"), covidx::QueryError);
}

TEST_CASE("CovIndex policies affect counts and manifest") {
  BuildOptions options;
  options.count_deletions = false;
  const CovIndex index = build_tiny(options);
  CHECK(index.total_covered("S2", "chr1") == 41);
  CHECK(index.covered_count("S2", "chr1", 200, 212) == 10);
  CHECK_FALSE(index.manifest().policy.count_deletions);
}

TEST_CASE("CovIndex groups and pooled units") {
  BuildOptions options;
  options.layout = Layout::Both;
  options.groups = {covidx::Group{"trio", {"S1", "S2"}}};
  const CovIndex index = build_tiny(options);
  REQUIRE(index.units().size() == 6);
  CHECK(index.units()[4].name == "trio");
  CHECK(index.units()[5].name == "pooled");
  CHECK(index.total_covered("trio", "chr1") == 103);
  CHECK(index.total_covered("pooled", "chr1") == 123);
  CHECK(index.covered_count("trio", "chr1", 0, 1000) == 103);
  CHECK(index.covered_count("pooled", "chr1", 0, 1000) == 123);
}

TEST_CASE("CovIndex pooled-only layout") {
  BuildOptions options;
  options.layout = Layout::Pooled;
  const CovIndex index = build_tiny(options);
  REQUIRE(index.units().size() == 1);
  CHECK(index.units()[0].name == "pooled");
  CHECK(index.total_covered("pooled", "chr1") == 123);
}

TEST_CASE("CovIndex scoped builds reject out-of-scope windows") {
  BuildOptions options;
  options.regions = {
      Region{"chr1", 100, 650},
      Region{"chr1", 900, 1000},
      Region{"chr2", 0, 500},
  };
  const CovIndex scoped = build_tiny(options);
  CHECK(scoped.manifest().index.partial_scope);
  REQUIRE(scoped.scope("chr1").size() == 2);
  CHECK(scoped.scope("chr1")[0].start == 100);
  CHECK(scoped.scope("chr1")[0].end == 650);
  CHECK(scoped.scope("chr1")[1].start == 900);
  CHECK(scoped.scope("chr1")[1].end == 1000);
  CHECK(scoped.total_covered("S1", "chr1") == 40);

  const CovIndex full = build_tiny();
  CHECK(scoped.covered_count("S1", "chr1", 100, 650) ==
        full.covered_count("S1", "chr1", 100, 650));
  CHECK(scoped.covered_count("S1", "chr1", 900, 1000) ==
        full.covered_count("S1", "chr1", 900, 1000));
  CHECK_THROWS_AS((void)scoped.covered_count("S1", "chr1", 0, 1000),
                  covidx::QueryError);
  CHECK_THROWS_AS((void)scoped.covered_count("S1", "chr1", 650, 900),
                  covidx::QueryError);
  CHECK_THROWS_AS((void)scoped.covered_count("S1", "chr1", 640, 660),
                  covidx::QueryError);
}

TEST_CASE("CovIndex save and load preserve answers") {
  tests::TempDir temp;
  CovIndex built = build_tiny();
  built.save(temp.path().string(), "tiny");
  const fs::path directory = temp.path() / "tiny.covidx";
  REQUIRE(fs::exists(directory / "manifest.json"));
  REQUIRE(fs::exists(directory / "data" / "u000000.bin"));
  REQUIRE(fs::exists(directory / "data" / "u000000.bin.crc32"));
  REQUIRE(fs::exists(directory / "contigs.tsv.crc32"));

  // Queries still work after save.
  CHECK(built.covered_count("S1", "chr1", 0, 1000) == 60);

  CovIndex loaded;
  loaded.load(directory.string());
  CHECK(loaded.covered_count("S1", "chr1", 0, 1000) == 60);
  CHECK_THROWS_AS((void)loaded.covered_count("pooled", "chr1", 0, 1),
                  covidx::QueryError);
  CHECK(loaded.manifest().index.unit_count == 4);
  CHECK(loaded.build_stats().records_accepted == 12);
  CHECK_THROWS_AS((void)loaded.save(temp.path().string(), "again"),
                  covidx::UsageError);

  // Heap mode returns the same answers.
  OpenOptions open;
  open.mmap = false;
  CovIndex heap;
  heap.load(directory.string(), open);
  CHECK(heap.covered_count("S2", "chr1", 200, 212) == 12);
}

TEST_CASE("CovIndex save overwrite and refusal rules") {
  tests::TempDir temp;
  CovIndex index = build_tiny();
  index.save(temp.path().string(), "tiny");
  CHECK_THROWS_AS((void)index.save(temp.path().string(), "tiny"),
                  covidx::OutputError);
  CHECK_NOTHROW(index.save(temp.path().string(), "tiny", true));
  CHECK_THROWS_AS((void)index.save(temp.path().string(), "bad/name", true),
                  covidx::UsageError);

  const fs::path target = temp.path() / "link.covidx";
  fs::create_symlink(temp.path() / "elsewhere", target);
  CHECK_THROWS_AS((void)index.save(temp.path().string(), "link", true),
                  covidx::OutputError);
}

TEST_CASE("CovIndex save refuses a corrupted work index") {
  tests::TempDir temp;
  BuildOptions options;
  options.work_dir = temp.path().string();
  options.keep_work = true;
  CovIndex index = build_tiny(options);
  const fs::path work = find_work_dir(temp.path());
  REQUIRE(!work.empty());
  {
    std::ofstream out(work / "data" / "u000000.bin",
                      std::ios::binary | std::ios::app);
    out.put('\0');
  }
  const fs::path out_dir = temp.path() / "out";
  CHECK_THROWS_AS((void)index.save(out_dir.string(), "tiny"),
                  covidx::OutputError);
  CHECK_FALSE(fs::exists(out_dir / "tiny.covidx"));
  bool temp_left = false;
  if (fs::exists(out_dir)) {
    for (const auto& entry : fs::directory_iterator(out_dir)) {
      if (entry.path().filename().string().find(".tmp.") != std::string::npos)
        temp_left = true;
    }
  }
  CHECK_FALSE(temp_left);
}

TEST_CASE("CovIndex work directory cleanup and keep_work") {
  tests::TempDir temp;
  const fs::path parent = temp.path() / "work";
  fs::create_directories(parent);
  BuildOptions options;
  options.work_dir = parent.string();
  {
    CovIndex index = build_tiny(options);
    const fs::path work = find_work_dir(parent);
    CHECK(!work.empty());
    // Finalization consumes the spools even though the work directory
    // remains alive for save().
    CHECK(fs::exists(work / "data"));
    CHECK_FALSE(fs::exists(work / "spools"));
    index.save(temp.path().string(), "tiny");
    OpenOptions open;
    index.load((temp.path() / "tiny.covidx").string(), open);
    CHECK(find_work_dir(parent).empty());
  }

  BuildOptions keep = options;
  keep.keep_work = true;
  {
    CovIndex index = build_tiny(keep);
    const fs::path work = find_work_dir(parent);
    CHECK(!work.empty());
    CHECK(fs::exists(work));
    CHECK(fs::exists(work / "spools"));
  }
  CHECK(!find_work_dir(parent).empty());
  fs::remove_all(parent);
}

TEST_CASE("CovIndex move semantics") {
  CovIndex index = build_tiny();
  CovIndex moved = std::move(index);
  CHECK_THROWS_AS((void)index.covered_count("S1", "chr1", 0, 10),
                  covidx::QueryError);
  CHECK(moved.covered_count("S1", "chr1", 0, 10) == 10);
}
