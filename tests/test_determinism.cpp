#include "doctest.h"

#include <cstdint>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include "covidx/detail/bam_reader.h"
#include "covidx/detail/merge.h"
#include "covidx/detail/spool.h"
#include "test_util.h"

namespace fs = std::filesystem;
using covidx::BuildOptions;
using covidx::BuildStats;
using covidx::Group;
using covidx::Layout;
using covidx::detail::BamHeaderInfo;
using covidx::detail::BuildWorkspace;
using covidx::detail::make_units;
using covidx::detail::read_bam_header;
using covidx::detail::UnitPlan;

namespace {

BamHeaderInfo tiny_header() {
  return read_bam_header(
      (fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam").string());
}

std::map<std::string, std::string> snapshot(const fs::path& root) {
  std::map<std::string, std::string> out;
  for (const auto& entry : fs::recursive_directory_iterator(root)) {
    if (!entry.is_regular_file())
      continue;
    out[fs::relative(entry.path(), root).generic_string()] =
        tests::read_file(entry.path());
  }
  return out;
}

BuildStats run_with_threads(const BamHeaderInfo& header,
                            const BuildOptions& base, unsigned threads,
                            BuildWorkspace& workspace,
                            const std::vector<UnitPlan>& units) {
  BuildOptions options = base;
  options.threads = threads;
  const std::vector<uint32_t> tids = {0, 1};
  const std::vector<std::vector<std::pair<uint64_t, uint64_t>>> intervals = {
      {{0, 1000}}, {{0, 500}}};
  return covidx::detail::build_run_spools(header, options, units, tids,
                                          intervals, workspace);
}

} // namespace

TEST_CASE("spool builds are byte-identical for 1/2/8 threads") {
  const BamHeaderInfo header = tiny_header();
  BuildOptions options;
  options.layout = Layout::Both;
  options.groups = {Group{"trio", {"S1", "S2"}}};
  const std::vector<UnitPlan> units = make_units(header, options);
  tests::TempDir parent;

  BuildWorkspace one(parent.path().string(), false);
  const BuildStats stats_one = run_with_threads(header, options, 1, one, units);
  const auto files_one = snapshot(one.root());

  for (unsigned threads : {2u, 8u}) {
    BuildWorkspace other(parent.path().string(), false);
    const BuildStats stats_other =
        run_with_threads(header, options, threads, other, units);
    CHECK(stats_other.records_examined == stats_one.records_examined);
    CHECK(stats_other.records_accepted == stats_one.records_accepted);
    CHECK(stats_other.runs == stats_one.runs);
    CHECK(snapshot(other.root()) == files_one);
  }
}
