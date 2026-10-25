#include "doctest.h"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <random>
#include <string>
#include <vector>

#include <htslib/sam.h>

#include "covidx/detail/bam_reader.h"
#include "covidx/detail/build.h"
#include "covidx/detail/merge.h"
#include "covidx/detail/run_index.h"
#include "covidx/error.h"
#include "test_util.h"

namespace fs = std::filesystem;
using covidx::AmbiguousBasePolicy;
using covidx::BuildOptions;
using covidx::Group;
using covidx::Layout;
using covidx::detail::BamHeaderInfo;
using covidx::detail::effective_exclude_flags;
using covidx::detail::make_units;
using covidx::detail::read_bam_header;
using covidx::detail::RunIndex;
using covidx::detail::UnitPlan;

namespace {

BamHeaderInfo tiny_header() {
  return read_bam_header(
      (fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam").string());
}

void add_synthetic_sample(BamHeaderInfo& header, const std::string& name) {
  if (header.sample_by_name.contains(name))
    return;
  const uint32_t id = static_cast<uint32_t>(header.samples.size());
  header.samples.push_back(covidx::detail::SampleRecord{name, {}});
  header.sample_by_name.emplace(name, id);
}

struct Bitmaps {
  std::vector<std::vector<bool>> samples;
  std::vector<bool> any;
};

Bitmaps naive_bitmaps(const BamHeaderInfo& header, const BuildOptions& options,
                      uint32_t tid) {
  const uint64_t length = header.contigs[tid].length;
  Bitmaps out;
  out.samples.assign(
      header.samples.size(),
      std::vector<bool>(static_cast<std::size_t>(length), false));
  out.any.assign(static_cast<std::size_t>(length), false);
  const uint16_t exclude = effective_exclude_flags(options);
  const bool demux = options.layout != Layout::Pooled;
  covidx::detail::BamRegionReader reader(header.path, static_cast<int32_t>(tid),
                                         {{0, length}});
  bam1_t* record = bam_init1();
  while (reader.next(*record)) {
    if ((record->core.flag & BAM_FUNMAP) != 0 || record->core.tid < 0 ||
        record->core.pos < 0) {
      continue;
    }
    if ((record->core.flag & exclude) != 0)
      continue;
    if (record->core.qual < options.min_mapq)
      continue;
    uint32_t sample = 0;
    if (demux)
      sample = covidx::detail::sample_id_for_read(*record, header,
                                                  options.unknown_sample);
    uint64_t ref = static_cast<uint64_t>(record->core.pos);
    uint64_t query = 0;
    const uint32_t* cigar = bam_get_cigar(record);
    const uint8_t* seq =
        record->core.l_qseq > 0 ? bam_get_seq(record) : nullptr;
    for (uint32_t i = 0; i < record->core.n_cigar; ++i) {
      const uint32_t len = cigar[i] >> 4;
      const uint32_t op = cigar[i] & 0x0F;
      if (op == 0 || op == 7 || op == 8) {
        for (uint32_t k = 0; k < len; ++k) {
          bool contributes = true;
          if (options.ambiguous_read_bases == AmbiguousBasePolicy::Exclude) {
            const uint8_t base = bam_seqi(seq, static_cast<int>(query + k));
            contributes = base == 1 || base == 2 || base == 4 || base == 8;
          }
          if (contributes) {
            out.any[static_cast<std::size_t>(ref + k)] = true;
            if (demux) {
              out.samples[sample][static_cast<std::size_t>(ref + k)] = true;
            }
          }
        }
        ref += len;
        query += len;
      } else if (op == 2) {
        if (options.count_deletions) {
          for (uint32_t k = 0; k < len; ++k) {
            out.any[static_cast<std::size_t>(ref + k)] = true;
            if (demux)
              out.samples[sample][static_cast<std::size_t>(ref + k)] = true;
          }
        }
        ref += len;
      } else if (op == 3) {
        ref += len;
      } else if (op == 1 || op == 4) {
        query += len;
      }
    }
  }
  bam_destroy1(record);
  return out;
}

std::vector<bool> unit_expectation(const UnitPlan& unit, const Bitmaps& bitmaps,
                                   std::size_t sample_count) {
  if (unit.kind == covidx::UnitInfo::Kind::Sample) {
    return bitmaps.samples.at(unit.sample_ids.at(0));
  }
  if (unit.kind == covidx::UnitInfo::Kind::Pooled && unit.sample_ids.empty()) {
    return bitmaps.any;
  }
  std::vector<bool> out(bitmaps.any.size(), false);
  for (uint32_t sample : unit.sample_ids) {
    REQUIRE(sample < sample_count);
    const std::vector<bool>& source = bitmaps.samples[sample];
    for (std::size_t i = 0; i < out.size(); ++i)
      out[i] = out[i] || source[i];
  }
  return out;
}

void compare_unit(const UnitPlan& unit,
                  const std::vector<covidx::detail::Run>& runs,
                  const std::vector<bool>& expected, uint64_t length) {
  (void)unit;
  RunIndex index;
  index.reset(length);
  index.set_runs(runs);
  uint64_t total = 0;
  for (bool value : expected)
    total += value ? 1 : 0;
  CHECK(index.total_covered() == total);
  for (uint64_t x = 0; x < length; ++x) {
    uint64_t running = 0;
    const uint64_t limit = std::min<uint64_t>(length, x + 60);
    for (uint64_t y = x + 1; y <= limit; ++y) {
      running += expected[static_cast<std::size_t>(y - 1)] ? 1 : 0;
      CHECK(index.covered_count(x, y) == running);
    }
  }
}

} // namespace

TEST_CASE("per-sample units match naive per-sample coverage") {
  const BamHeaderInfo header = tiny_header();
  BuildOptions options;
  const std::vector<UnitPlan> units = make_units(header, options);
  REQUIRE(units.size() == 4);
  CHECK(units[0].name == "S1");
  CHECK(units[1].name == "S2");
  CHECK(units[2].name == "S3");
  CHECK(units[3].name == "S4");
  CHECK(units[2].kind == covidx::UnitInfo::Kind::Sample);

  std::vector<std::vector<covidx::detail::Run>> runs;
  covidx::detail::build_contig_all_units(header, options, units, 0,
                                         {{0, header.contigs[0].length}}, runs);
  const Bitmaps bitmaps = naive_bitmaps(header, options, 0);
  for (const UnitPlan& unit : units) {
    compare_unit(unit, runs.at(unit.id),
                 unit_expectation(unit, bitmaps, header.samples.size()),
                 header.contigs[0].length);
  }
}

TEST_CASE("pooled and fixed groups are unions, not sums") {
  const BamHeaderInfo header = tiny_header();
  BuildOptions options;
  options.layout = Layout::Both;
  options.groups = {Group{"trio", {"S1", "S2"}},
                    Group{"everyone", {"S1", "S2", "S3", "S4"}}};
  const std::vector<UnitPlan> units = make_units(header, options);
  REQUIRE(units.size() == 7);
  CHECK(units[4].name == "everyone");
  CHECK(units[4].kind == covidx::UnitInfo::Kind::Group);
  CHECK(units[5].name == "trio");
  CHECK(units[6].name == "pooled");
  CHECK(units[6].kind == covidx::UnitInfo::Kind::Pooled);
  CHECK(units[6].members.size() == 4);

  std::vector<std::vector<covidx::detail::Run>> runs;
  covidx::detail::build_contig_all_units(header, options, units, 0,
                                         {{0, header.contigs[0].length}}, runs);
  const Bitmaps bitmaps = naive_bitmaps(header, options, 0);
  for (const UnitPlan& unit : units) {
    compare_unit(unit, runs.at(unit.id),
                 unit_expectation(unit, bitmaps, header.samples.size()),
                 header.contigs[0].length);
  }

  // The pooled union must be exactly the union of all sample runs.
  std::vector<const std::vector<covidx::detail::Run>*> sample_lists;
  for (std::size_t i = 0; i < header.samples.size(); ++i) {
    sample_lists.push_back(&runs.at(i));
  }
  const auto merged = covidx::detail::merge_run_lists(sample_lists);
  CHECK(merged.size() == runs.at(units.back().id).size());
}

TEST_CASE(
    "pooled-only layout ignores RG and matches a full single-unit build") {
  const BamHeaderInfo header = tiny_header();
  BuildOptions options;
  options.layout = Layout::Pooled;
  const std::vector<UnitPlan> units = make_units(header, options);
  REQUIRE(units.size() == 1);
  CHECK(units[0].name == "pooled");
  CHECK(units[0].direct_from_reads);

  std::vector<std::vector<covidx::detail::Run>> runs;
  covidx::detail::build_contig_all_units(header, options, units, 0,
                                         {{0, header.contigs[0].length}}, runs);
  RunIndex single;
  covidx::detail::build_single_unit(header, options, 0, single);
  RunIndex pooled;
  pooled.reset(header.contigs[0].length);
  pooled.set_runs(runs.at(0));
  CHECK(pooled.total_covered() == single.total_covered());
  for (uint64_t x = 0; x < 1000; x += 17) {
    for (uint64_t y = x + 1; y < 1000; y += 53) {
      CHECK(pooled.covered_count(x, y) == single.covered_count(x, y));
    }
  }
}

TEST_CASE("unknown-sample reads join a synthetic sample") {
  tests::TempDir temp;
  const std::string header_text = "@HD\tVN:1.6\tSO:coordinate\n"
                                  "@SQ\tSN:chr1\tLN:500\n"
                                  "@RG\tID:rg1\tSM:S1\n"
                                  "@RG\tID:rg2\tSM:S2\n";
  std::vector<tests::ReadSpec> reads;
  reads.push_back(
      {"known", 0, 0, 10, 60, {{'M', 10}}, tests::repeat_sequence(10), "rg1"});
  reads.push_back(
      {"no_rg", 0, 0, 100, 60, {{'M', 10}}, tests::repeat_sequence(10), ""});
  reads.push_back({"bad_rg",
                   0,
                   0,
                   200,
                   60,
                   {{'M', 10}},
                   tests::repeat_sequence(10),
                   "rgX"});
  const fs::path bam = temp.path() / "unknown.bam";
  tests::make_bam(bam, header_text, reads);

  BamHeaderInfo header = read_bam_header(bam.string());
  add_synthetic_sample(header, "UNK");
  BuildOptions options;
  options.unknown_sample = "UNK";
  const std::vector<UnitPlan> units = make_units(header, options);
  REQUIRE(units.size() == 3);
  CHECK(units[2].name == "UNK");

  std::vector<std::vector<covidx::detail::Run>> runs;
  covidx::detail::build_contig_all_units(header, options, units, 0, {{0, 500}},
                                         runs);
  const Bitmaps bitmaps = naive_bitmaps(header, options, 0);
  for (const UnitPlan& unit : units) {
    compare_unit(unit, runs.at(unit.id),
                 unit_expectation(unit, bitmaps, header.samples.size()), 500);
  }
  RunIndex unknown;
  unknown.reset(500);
  unknown.set_runs(runs.at(2));
  CHECK(unknown.covered_count(100, 110) == 10);
  CHECK(unknown.covered_count(200, 210) == 10);
  CHECK(unknown.covered_count(0, 10) == 0);
}

TEST_CASE("unit planning validation matrix") {
  const BamHeaderInfo header = tiny_header();

  BuildOptions pooled_group;
  pooled_group.layout = Layout::Pooled;
  pooled_group.groups = {Group{"g", {"S1"}}};
  CHECK_THROWS_AS((void)make_units(header, pooled_group), covidx::UsageError);

  BuildOptions unknown_member;
  unknown_member.groups = {Group{"g", {"S1", "nope"}}};
  CHECK_THROWS_AS((void)make_units(header, unknown_member), covidx::UsageError);

  BuildOptions duplicate_member;
  duplicate_member.groups = {Group{"g", {"S1", "S1"}}};
  CHECK_THROWS_AS((void)make_units(header, duplicate_member),
                  covidx::UsageError);

  BuildOptions collision;
  collision.groups = {Group{"S1", {"S2"}}};
  CHECK_THROWS_AS((void)make_units(header, collision), covidx::UsageError);

  BuildOptions reserved;
  reserved.groups = {Group{"pooled", {"S1"}}};
  CHECK_THROWS_AS((void)make_units(header, reserved), covidx::UsageError);

  BuildOptions duplicate_group;
  duplicate_group.groups = {Group{"g", {"S1"}}, Group{"g", {"S2"}}};
  CHECK_THROWS_AS((void)make_units(header, duplicate_group),
                  covidx::UsageError);

  BuildOptions empty_group;
  empty_group.groups = {Group{"g", {}}};
  CHECK_THROWS_AS((void)make_units(header, empty_group), covidx::UsageError);
}
