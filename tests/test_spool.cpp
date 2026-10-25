#include "doctest.h"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "covidx/detail/bam_reader.h"
#include "covidx/detail/build.h"
#include "covidx/detail/merge.h"
#include "covidx/detail/run_index.h"
#include "covidx/detail/spool.h"
#include "covidx/error.h"
#include "test_util.h"

namespace fs = std::filesystem;
using covidx::BuildOptions;
using covidx::BuildStats;
using covidx::Group;
using covidx::Layout;
using covidx::detail::BamHeaderInfo;
using covidx::detail::BuildWorkspace;
using covidx::detail::make_units;
using covidx::detail::merge_run_spools;
using covidx::detail::read_bam_header;
using covidx::detail::Run;
using covidx::detail::RunSpoolReader;
using covidx::detail::RunSpoolWriter;
using covidx::detail::spool_path_for_unit;
using covidx::detail::UnitPlan;

namespace {

BamHeaderInfo tiny_header() {
  return read_bam_header(
      (fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam").string());
}

void check_runs(const std::vector<Run>& got, const std::vector<Run>& want) {
  REQUIRE(got.size() == want.size());
  for (std::size_t i = 0; i < want.size(); ++i) {
    CHECK(got[i].start == want[i].start);
    CHECK(got[i].end == want[i].end);
  }
}

std::map<std::string, std::string> snapshot(const fs::path& root) {
  std::map<std::string, std::string> out;
  for (const auto& entry : fs::recursive_directory_iterator(root)) {
    if (!entry.is_regular_file())
      continue;
    const std::string relative =
        fs::relative(entry.path(), root).generic_string();
    out[relative] = tests::read_file(entry.path());
  }
  return out;
}

} // namespace

TEST_CASE("spool round-trip and validation") {
  tests::TempDir temp;
  const fs::path path = temp.path() / "runs.bin";
  {
    RunSpoolWriter writer(path.string());
    writer.write(Run{0, 5});
    writer.write(Run{10, 20});
    writer.close();
    CHECK(writer.records() == 2);
  }
  RunSpoolReader reader(path.string());
  check_runs(reader.read_all(), {Run{0, 5}, Run{10, 20}});
  CHECK(fs::file_size(path) == 32);

  RunSpoolReader missing((temp.path() / "missing.bin").string());
  CHECK(missing.read_all().empty());

  const fs::path truncated = temp.path() / "truncated.bin";
  tests::write_file(truncated, std::string(8, '\0'));
  RunSpoolReader bad(truncated.string());
  Run run;
  CHECK_THROWS_AS((void)bad.next(run), covidx::FormatError);

  const fs::path unsorted = temp.path() / "unsorted.bin";
  {
    RunSpoolWriter writer(unsorted.string());
    writer.write(Run{10, 20});
    writer.write(Run{5, 8});
    writer.close();
  }
  RunSpoolReader bad_order(unsorted.string());
  CHECK_THROWS_AS((void)bad_order.read_all(), covidx::FormatError);

  const fs::path partial = temp.path() / "partial.bin";
  {
    RunSpoolWriter writer(partial.string());
    for (uint64_t i = 0; i < 4096; ++i)
      writer.write(Run{i * 10, i * 10 + 5});
    writer.close();
  }
  {
    std::ofstream out(partial, std::ios::binary | std::ios::app);
    out << std::string(8, '\0');
  }
  RunSpoolReader bad_partial(partial.string());
  CHECK_THROWS_AS((void)bad_partial.read_all(), covidx::FormatError);
}

TEST_CASE("feature spool round-trip and validation") {
  tests::TempDir temp;
  const fs::path path = temp.path() / "features.bin";
  {
    covidx::detail::FeatureSpoolWriter writer(path.string());
    writer.write(covidx::detail::FeatureSpan{0, 5});
    writer.write(covidx::detail::FeatureSpan{3, 10});
    writer.write(covidx::detail::FeatureSpan{20, 30});
    writer.close();
    CHECK(writer.records() == 3);
  }
  const std::vector<covidx::detail::FeatureSpan> spans =
      covidx::detail::FeatureSpoolReader(path.string()).read_all();
  REQUIRE(spans.size() == 3);
  CHECK(spans[1].start == 3);
  CHECK(spans[1].end == 10);
  CHECK(fs::file_size(path) == 24);

  covidx::detail::FeatureSpoolReader missing(
      (temp.path() / "missing.bin").string());
  CHECK(missing.read_all().empty());

  const fs::path truncated = temp.path() / "truncated.bin";
  tests::write_file(truncated, std::string(4, '\0'));
  CHECK_THROWS_AS(
      (void)covidx::detail::FeatureSpoolReader(truncated.string()).read_all(),
      covidx::FormatError);

  const fs::path unsorted = temp.path() / "unsorted.bin";
  {
    covidx::detail::FeatureSpoolWriter writer(unsorted.string());
    writer.write(covidx::detail::FeatureSpan{10, 20});
    writer.write(covidx::detail::FeatureSpan{5, 8});
    writer.close();
  }
  CHECK_THROWS_AS(
      (void)covidx::detail::FeatureSpoolReader(unsorted.string()).read_all(),
      covidx::FormatError);

  const fs::path partial = temp.path() / "feature_partial.bin";
  {
    covidx::detail::FeatureSpoolWriter writer(partial.string());
    for (uint64_t i = 0; i < 8192; ++i)
      writer.write(covidx::detail::FeatureSpan{
          static_cast<uint32_t>(i * 10), static_cast<uint32_t>(i * 10 + 5)});
    writer.close();
  }
  {
    std::ofstream out(partial, std::ios::binary | std::ios::app);
    out << std::string(4, '\0');
  }
  CHECK_THROWS_AS(
      (void)covidx::detail::FeatureSpoolReader(partial.string()).read_all(),
      covidx::FormatError);
}

TEST_CASE("merge_run_spools is an exact union") {
  tests::TempDir temp;
  const fs::path a = temp.path() / "a.bin";
  const fs::path b = temp.path() / "b.bin";
  const fs::path out = temp.path() / "out.bin";
  {
    RunSpoolWriter writer(a.string());
    writer.write(Run{0, 5});
    writer.write(Run{10, 20});
    writer.close();
  }
  {
    RunSpoolWriter writer(b.string());
    writer.write(Run{3, 12});
    writer.write(Run{25, 30});
    writer.close();
  }
  {
    std::vector<std::unique_ptr<RunSpoolReader>> readers;
    readers.push_back(std::make_unique<RunSpoolReader>(a.string()));
    readers.push_back(std::make_unique<RunSpoolReader>(b.string()));
    RunSpoolWriter writer(out.string());
    merge_run_spools(readers, writer);
    writer.close();
  }
  RunSpoolReader reader(out.string());
  check_runs(reader.read_all(), {Run{0, 20}, Run{25, 30}});
}

TEST_CASE("workspace owns and optionally keeps its directory") {
  tests::TempDir parent;
  fs::path root;
  {
    BuildWorkspace workspace(parent.path().string(), false);
    root = workspace.root();
    CHECK(fs::exists(root));
    CHECK(fs::exists(workspace.data_dir()));
    CHECK(fs::exists(workspace.spool_dir()));
    CHECK(workspace.sample_spool(3, 7).filename() == "c000007.bin");
    CHECK(workspace.unit_spool(3, 7).parent_path().filename() == "u000003");
  }
  CHECK_FALSE(fs::exists(root));

  {
    BuildWorkspace workspace(parent.path().string(), true);
    root = workspace.root();
    CHECK(fs::exists(root));
  }
  CHECK(fs::exists(root));
  fs::remove_all(root);
}

TEST_CASE("build_run_spools matches in-memory per-sample builds") {
  const BamHeaderInfo header = tiny_header();
  BuildOptions options;
  options.threads = 2;
  const std::vector<UnitPlan> units = make_units(header, options);
  tests::TempDir parent;
  BuildWorkspace workspace(parent.path().string(), false);
  const std::vector<uint32_t> tids = {0, 1};
  const std::vector<std::vector<std::pair<uint64_t, uint64_t>>> intervals = {
      {{0, 1000}}, {{0, 500}}};
  const BuildStats stats = covidx::detail::build_run_spools(
      header, options, units, tids, intervals, workspace);
  CHECK(stats.units == 4);
  CHECK(stats.contigs == 2);
  CHECK(stats.records_examined == tests::tiny_reads().size());

  for (std::size_t contig = 0; contig < tids.size(); ++contig) {
    std::vector<std::vector<Run>> in_memory;
    covidx::detail::build_contig_all_units(header, options, units, tids[contig],
                                           intervals[contig], in_memory);
    for (const UnitPlan& unit : units) {
      RunSpoolReader reader(
          spool_path_for_unit(workspace, unit, static_cast<uint32_t>(contig))
              .string());
      check_runs(reader.read_all(), in_memory.at(unit.id));
    }
  }
}

TEST_CASE("build_run_spools handles groups, pooled, and pooled-only layouts") {
  const BamHeaderInfo header = tiny_header();
  tests::TempDir parent;
  const std::vector<uint32_t> tids = {0, 1};
  const std::vector<std::vector<std::pair<uint64_t, uint64_t>>> intervals = {
      {{0, 1000}}, {{0, 500}}};

  for (Layout layout : {Layout::Both, Layout::Pooled}) {
    BuildOptions options;
    options.layout = layout;
    options.threads = 2;
    if (layout == Layout::Both) {
      options.groups = {Group{"trio", {"S1", "S2"}}};
    }
    const std::vector<UnitPlan> units = make_units(header, options);
    BuildWorkspace workspace(parent.path().string(), false);
    (void)covidx::detail::build_run_spools(header, options, units, tids,
                                           intervals, workspace);
    for (std::size_t contig = 0; contig < tids.size(); ++contig) {
      std::vector<std::vector<Run>> in_memory;
      covidx::detail::build_contig_all_units(
          header, options, units, tids[contig], intervals[contig], in_memory);
      for (const UnitPlan& unit : units) {
        RunSpoolReader reader(
            spool_path_for_unit(workspace, unit, static_cast<uint32_t>(contig))
                .string());
        check_runs(reader.read_all(), in_memory.at(unit.id));
      }
    }
  }
}

TEST_CASE("sample batching produces identical spools and stats") {
  const BamHeaderInfo header = tiny_header();
  BuildOptions options;
  options.threads = 1;
  const std::vector<UnitPlan> units = make_units(header, options);
  const std::vector<uint32_t> tids = {0, 1};
  const std::vector<std::vector<std::pair<uint64_t, uint64_t>>> intervals = {
      {{0, 1000}}, {{0, 500}}};
  tests::TempDir parent;
  BuildWorkspace one(parent.path().string(), false);
  const BuildStats stats_one = covidx::detail::build_run_spools(
      header, options, units, tids, intervals, one, 1);
  BuildWorkspace many(parent.path().string(), false);
  const BuildStats stats_many = covidx::detail::build_run_spools(
      header, options, units, tids, intervals, many, 256);
  CHECK(stats_one.records_examined == stats_many.records_examined);
  CHECK(stats_one.records_accepted == stats_many.records_accepted);
  CHECK(stats_one.runs == stats_many.runs);
  CHECK(snapshot(one.root()) == snapshot(many.root()));
}

TEST_CASE("spool write failures abort the build and clean the workspace") {
  const BamHeaderInfo header = tiny_header();
  BuildOptions options;
  options.threads = 1;
  const std::vector<UnitPlan> units = make_units(header, options);
  tests::TempDir parent;
  fs::path root;
  {
    BuildWorkspace workspace(parent.path().string(), false);
    root = workspace.root();
    tests::write_file(workspace.spool_dir() / "s000000",
                      "block the sample directory");
    CHECK_THROWS_AS((void)covidx::detail::build_run_spools(
                        header, options, units, {0}, {{{0, 1000}}}, workspace),
                    covidx::OutputError);
  }
  CHECK_FALSE(fs::exists(root));
}
