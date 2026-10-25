#include "doctest.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <map>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <htslib/sam.h>

#include "covidx/cov_index.h"
#include "covidx/detail/bam_reader.h"
#include "covidx/detail/build.h"
#include "covidx/error.h"
#include "covidx/manifest.h"
#include "test_util.h"

namespace fs = std::filesystem;
using covidx::BuildOptions;
using covidx::CovIndex;
using covidx::Layout;
using covidx::QuerySpan;
using covidx::Region;
using covidx::detail::BamHeaderInfo;
using covidx::detail::BamRegionReader;
using covidx::detail::effective_exclude_flags;
using covidx::detail::read_bam_header;

namespace {

using Span = std::pair<uint64_t, uint64_t>;

fs::path tiny_bam() {
  return fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam";
}

CovIndex build_tiny(const BuildOptions& options = {}) {
  CovIndex index;
  index.build(tiny_bam().string(), options);
  return index;
}

uint32_t read_sample(const bam1_t& read, const BamHeaderInfo& header,
                     const BuildOptions& options) {
  if (options.layout == Layout::Pooled)
    return 0;
  return covidx::detail::sample_id_for_read(read, header,
                                            options.unknown_sample);
}

std::vector<Span> naive_spans(const BamHeaderInfo& header,
                              const BuildOptions& options, uint32_t tid,
                              int32_t sample_filter = -1) {
  const uint64_t length = header.contigs.at(tid).length;
  std::vector<Span> spans;
  BamRegionReader reader(header.path, static_cast<int32_t>(tid), {{0, length}});
  bam1_t* record = bam_init1();
  const uint16_t exclude = effective_exclude_flags(options);
  while (reader.next(*record)) {
    if ((record->core.flag & BAM_FUNMAP) != 0 || record->core.tid < 0 ||
        record->core.pos < 0)
      continue;
    if ((record->core.flag & exclude) != 0)
      continue;
    if (record->core.qual < options.min_mapq)
      continue;
    if (sample_filter >= 0 && read_sample(*record, header, options) !=
                                  static_cast<uint32_t>(sample_filter))
      continue;
    uint64_t span = 0;
    const uint32_t* cigar = bam_get_cigar(record);
    for (uint32_t i = 0; i < record->core.n_cigar; ++i) {
      const uint32_t op = cigar[i] & 0x0F;
      if (op == 0 || op == 2 || op == 3 || op == 7 || op == 8)
        span += cigar[i] >> 4;
    }
    if (span == 0)
      continue;
    const uint64_t start = static_cast<uint64_t>(record->core.pos);
    spans.emplace_back(start, start + span);
  }
  bam_destroy1(record);
  return spans;
}

class SourceDateEpoch {
public:
  SourceDateEpoch() {
    ::setenv("SOURCE_DATE_EPOCH", "0", 1);
  }
  ~SourceDateEpoch() {
    ::unsetenv("SOURCE_DATE_EPOCH");
  }
  SourceDateEpoch(const SourceDateEpoch&) = delete;
  SourceDateEpoch& operator=(const SourceDateEpoch&) = delete;
};

std::map<std::string, std::string> snapshot(const fs::path& root) {
  std::map<std::string, std::string> files;
  for (const auto& entry : fs::recursive_directory_iterator(root)) {
    if (!entry.is_regular_file())
      continue;
    files[fs::relative(entry.path(), root).generic_string()] =
        tests::read_file(entry.path());
  }
  return files;
}

uint64_t overlap_count(const std::vector<Span>& spans, uint64_t x, uint64_t y) {
  uint64_t count = 0;
  for (const Span& span : spans) {
    if (span.first < y && span.second > x)
      ++count;
  }
  return count;
}

void check_index_against(const CovIndex& index, const std::string& unit,
                         const BamHeaderInfo& header,
                         const BuildOptions& options, uint32_t tid,
                         int32_t sample_filter, uint32_t seed) {
  const std::vector<Span> spans =
      naive_spans(header, options, tid, sample_filter);
  const uint64_t length = header.contigs.at(tid).length;
  const std::string& contig = header.contigs.at(tid).name;
  CHECK(index.feature_count(unit, contig, 0, length) == spans.size());
  std::mt19937 rng(seed);
  for (int i = 0; i < 500; ++i) {
    const uint64_t x = rng() % length;
    const uint64_t y = x + 1 + rng() % (length - x);
    CHECK(index.feature_count(unit, contig, x, y) ==
          overlap_count(spans, x, y));
  }
}

} // namespace

TEST_CASE("feature counts match naive alignment spans") {
  const CovIndex index = build_tiny();
  CHECK(index.manifest().format_version == "1.1.0");
  REQUIRE(index.units().size() == 4);
  CHECK(index.units()[0].total_features == 6);
  CHECK(index.units()[1].total_features == 4);
  CHECK(index.units()[2].total_features == 0);
  CHECK(index.units()[3].total_features == 2);
  CHECK(index.build_stats().features == 12);

  const BamHeaderInfo header = read_bam_header(tiny_bam().string());
  BuildOptions options;
  for (std::size_t sample = 0; sample < header.samples.size(); ++sample) {
    const std::string& name = header.samples[sample].name;
    check_index_against(index, name, header, options, 0,
                        static_cast<int32_t>(sample),
                        static_cast<uint32_t>(100 + sample));
  }

  CHECK(index.feature_count("S1", "chr1", 400, 800) == 3);
  CHECK(index.feature_count("S2", "chr1", 200, 300) == 1);
  CHECK(index.feature_count("S2", "chr1", 200, 212) == 1);
  CHECK(index.feature_count("S3", "chr1", 0, 1000) == 0);
  CHECK(index.feature_count("S1", "chr2", 0, 500) == 0);
  CHECK(index.feature_count("S1", "chr1", 0, 10) == 1);
}

TEST_CASE("group and pooled feature counts sum member samples") {
  BuildOptions options;
  options.layout = Layout::Both;
  options.groups = {covidx::Group{"trio", {"S1", "S2"}}};
  const CovIndex index = build_tiny(options);
  REQUIRE(index.units().size() == 6);
  CHECK(index.unit("trio")->total_features == 10);
  CHECK(index.unit("pooled")->total_features == 12);
  CHECK(index.build_stats().features == 12);

  const BamHeaderInfo header = read_bam_header(tiny_bam().string());
  const uint64_t length = header.contigs[0].length;
  std::vector<Span> trio;
  std::vector<Span> pooled;
  for (std::size_t sample = 0; sample < header.samples.size(); ++sample) {
    const std::vector<Span> spans =
        naive_spans(header, options, 0, static_cast<int32_t>(sample));
    pooled.insert(pooled.end(), spans.begin(), spans.end());
    if (sample < 2)
      trio.insert(trio.end(), spans.begin(), spans.end());
  }
  std::mt19937 rng(901);
  for (int i = 0; i < 500; ++i) {
    const uint64_t x = rng() % length;
    const uint64_t y = x + 1 + rng() % (length - x);
    CHECK(index.feature_count("trio", "chr1", x, y) ==
          overlap_count(trio, x, y));
    CHECK(index.feature_count("pooled", "chr1", x, y) ==
          overlap_count(pooled, x, y));
  }
}

TEST_CASE("batch feature counts match single queries") {
  BuildOptions options;
  options.layout = Layout::Both;
  options.groups = {covidx::Group{"trio", {"S1", "S2"}}};
  const CovIndex index = build_tiny(options);
  const uint64_t length = 1000;
  std::mt19937 rng(711);
  std::vector<QuerySpan> spans;
  for (int i = 0; i < 2000; ++i) {
    const uint64_t x = rng() % length;
    const uint64_t y = x + 1 + rng() % (length - x);
    spans.push_back(QuerySpan{x, y});
  }
  for (const char* unit : {"S1", "S2", "trio", "pooled"}) {
    const std::vector<uint64_t> batch =
        index.feature_counts(unit, "chr1", spans);
    REQUIRE(batch.size() == spans.size());
    for (std::size_t i = 0; i < spans.size(); ++i) {
      CHECK(batch[i] ==
            index.feature_count(unit, "chr1", spans[i].start, spans[i].end));
    }
  }
  const std::vector<QuerySpan> no_spans;
  const std::vector<uint64_t> empty =
      index.feature_counts("S1", "chr1", no_spans);
  CHECK(empty.empty());
  const std::vector<QuerySpan> chr2{QuerySpan{0, 100}, QuerySpan{100, 200}};
  const std::vector<uint64_t> no_features =
      index.feature_counts("S1", "chr2", chr2);
  CHECK(no_features == std::vector<uint64_t>{0, 0});
  CHECK_THROWS_AS((void)index.feature_counts("nope", "chr1", spans),
                  covidx::QueryError);
  CHECK_THROWS_AS((void)index.feature_counts("S1", "nope", spans),
                  covidx::QueryError);
  const std::vector<QuerySpan> reversed{QuerySpan{5, 5}};
  CHECK_THROWS_AS((void)index.feature_counts("S1", "chr1", reversed),
                  covidx::QueryError);
  const std::vector<QuerySpan> beyond{QuerySpan{0, 1001}};
  CHECK_THROWS_AS((void)index.feature_counts("S1", "chr1", beyond),
                  covidx::QueryError);
}

TEST_CASE("pooled-only layout stores features directly") {
  BuildOptions options;
  options.layout = Layout::Pooled;
  const CovIndex index = build_tiny(options);
  REQUIRE(index.units().size() == 1);
  CHECK(index.unit("pooled")->total_features == 12);
  CHECK(index.feature_count("pooled", "chr1", 0, 1000) == 12);
  CHECK(index.feature_count("pooled", "chr1", 0, 100) == 1);
}

TEST_CASE("feature spans ignore deletion and ambiguous-base coverage policy") {
  const CovIndex defaults = build_tiny();
  BuildOptions no_deletions;
  no_deletions.count_deletions = false;
  const CovIndex gaps = build_tiny(no_deletions);
  CHECK(defaults.covered_count("S2", "chr1", 200, 212) == 12);
  CHECK(gaps.covered_count("S2", "chr1", 200, 212) == 10);
  CHECK(defaults.feature_count("S2", "chr1", 200, 212) == 1);
  CHECK(gaps.feature_count("S2", "chr1", 200, 212) == 1);

  tests::TempDir temp;
  const std::string header_text = "@HD\tVN:1.6\tSO:coordinate\n"
                                  "@SQ\tSN:chr1\tLN:200\n"
                                  "@RG\tID:rg1\tSM:S1\n";
  std::vector<tests::ReadSpec> reads;
  reads.push_back({"ambig", 0, 0, 0, 60, {{'M', 10}}, "ACGTNNACGT", "rg1"});
  reads.push_back(
      {"all_ambig", 0, 0, 50, 60, {{'M', 10}}, "NNNNNNNNNN", "rg1"});
  reads.push_back(
      {"plain", 0, 0, 100, 60, {{'M', 10}}, tests::repeat_sequence(10), "rg1"});
  const fs::path bam = temp.path() / "ambig.bam";
  tests::make_bam(bam, header_text, reads);

  BuildOptions cover;
  CovIndex cover_index;
  cover_index.build(bam.string(), cover);
  BuildOptions exclude = cover;
  exclude.ambiguous_read_bases = covidx::AmbiguousBasePolicy::Exclude;
  CovIndex exclude_index;
  exclude_index.build(bam.string(), exclude);
  CHECK(cover_index.feature_count("S1", "chr1", 0, 10) == 1);
  CHECK(exclude_index.feature_count("S1", "chr1", 0, 10) == 1);
  CHECK(cover_index.feature_count("S1", "chr1", 50, 60) == 1);
  CHECK(exclude_index.feature_count("S1", "chr1", 50, 60) == 1);
  CHECK(cover_index.feature_count("S1", "chr1", 0, 200) == 3);
  CHECK(exclude_index.feature_count("S1", "chr1", 0, 200) == 3);
  CHECK(cover_index.covered_count("S1", "chr1", 0, 10) == 10);
  CHECK(exclude_index.covered_count("S1", "chr1", 0, 10) == 8);
  CHECK(cover_index.covered_count("S1", "chr1", 50, 60) == 10);
  CHECK(exclude_index.covered_count("S1", "chr1", 50, 60) == 0);
}

TEST_CASE("feature spans include N and skip zero-reference-span reads") {
  tests::TempDir temp;
  const std::string header_text = "@HD\tVN:1.6\tSO:coordinate\n"
                                  "@SQ\tSN:chr1\tLN:500\n"
                                  "@RG\tID:rg1\tSM:S1\n";
  std::vector<tests::ReadSpec> reads;
  reads.push_back(
      {"normal", 0, 0, 0, 60, {{'M', 10}}, tests::repeat_sequence(10), "rg1"});
  reads.push_back({"all_n", 0, 0, 100, 60, {{'N', 10}}, "", "rg1"});
  reads.push_back({"zero", 0, 0, 200, 60, {{'S', 10}}, "ACGTACGTAC", "rg1"});
  const fs::path bam = temp.path() / "spans.bam";
  tests::make_bam(bam, header_text, reads);

  BuildOptions options;
  CovIndex index;
  index.build(bam.string(), options);
  CHECK(index.unit("S1")->total_features == 2);
  CHECK(index.feature_count("S1", "chr1", 0, 10) == 1);
  CHECK(index.feature_count("S1", "chr1", 100, 110) == 1);
  CHECK(index.feature_count("S1", "chr1", 200, 210) == 0);
  CHECK(index.feature_count("S1", "chr1", 0, 500) == 2);
}

TEST_CASE("feature counts span multiple contigs") {
  tests::TempDir temp;
  const std::string header_text = "@HD\tVN:1.6\tSO:coordinate\n"
                                  "@SQ\tSN:chr1\tLN:100\n"
                                  "@SQ\tSN:chr2\tLN:200\n"
                                  "@RG\tID:rg1\tSM:S1\n";
  std::vector<tests::ReadSpec> reads;
  reads.push_back(
      {"r1", 0, 0, 0, 60, {{'M', 10}}, tests::repeat_sequence(10), "rg1"});
  reads.push_back({"r2",
                   0,
                   1,
                   50,
                   60,
                   {{'M', 5}, {'N', 5}, {'M', 5}},
                   tests::repeat_sequence(10),
                   "rg1"});
  const fs::path bam = temp.path() / "two.bam";
  tests::make_bam(bam, header_text, reads);

  BuildOptions options;
  CovIndex index;
  index.build(bam.string(), options);
  CHECK(index.unit("S1")->total_features == 2);
  CHECK(index.feature_count("S1", "chr1", 0, 10) == 1);
  CHECK(index.feature_count("S1", "chr1", 10, 100) == 0);
  CHECK(index.feature_count("S1", "chr2", 50, 65) == 1);
  CHECK(index.feature_count("S1", "chr2", 0, 50) == 0);
}

TEST_CASE("scoped feature counts keep full alignment spans") {
  const CovIndex full = build_tiny();

  BuildOptions scoped_options;
  scoped_options.regions = {Region{"chr1", 205, 206}};
  const CovIndex scoped = build_tiny(scoped_options);
  CHECK(scoped.unit("S2")->total_features == 1);
  CHECK(scoped.feature_count("S2", "chr1", 205, 206) == 1);
  CHECK(scoped.feature_count("S2", "chr1", 205, 206) ==
        full.feature_count("S2", "chr1", 205, 206));

  BuildOptions scope_options;
  scope_options.regions = {
      Region{"chr1", 100, 350},
      Region{"chr1", 500, 650},
      Region{"chr1", 900, 1000},
  };
  const CovIndex scope = build_tiny(scope_options);
  std::mt19937 rng(7331);
  for (const Region& interval : scope_options.regions) {
    for (int i = 0; i < 200; ++i) {
      const uint64_t x =
          interval.start + rng() % (interval.end - interval.start);
      const uint64_t y = x + 1 + rng() % (interval.end - x);
      for (const char* unit : {"S1", "S2", "S3", "S4"}) {
        CHECK(scope.feature_count(unit, "chr1", x, y) ==
              full.feature_count(unit, "chr1", x, y));
      }
    }
  }
}

TEST_CASE("feature counts survive save, load, mmap, and heap") {
  tests::TempDir temp;
  BuildOptions options;
  {
    CovIndex index = build_tiny(options);
    index.save(temp.path().string(), "features");
  }
  const fs::path directory = temp.path() / "features.covidx";

  covidx::OpenOptions mmap_options;
  mmap_options.mmap = true;
  covidx::OpenOptions heap_options;
  heap_options.mmap = false;
  CovIndex mmap_index;
  mmap_index.load(directory.string(), mmap_options);
  CovIndex heap_index;
  heap_index.load(directory.string(), heap_options);

  std::mt19937 rng(4242);
  for (int i = 0; i < 2000; ++i) {
    const std::string unit = i % 2 == 0 ? "S1" : "S2";
    const uint64_t x = rng() % 1000;
    const uint64_t y = x + 1 + rng() % (1000 - x);
    CHECK(mmap_index.feature_count(unit, "chr1", x, y) ==
          heap_index.feature_count(unit, "chr1", x, y));
    CHECK(mmap_index.feature_count(unit, "chr1", x, y) ==
          heap_index.feature_count(unit, "chr1", x, y));
  }
  CHECK(mmap_index.unit("S1")->total_features == 6);
  CHECK(heap_index.unit("pooled") == nullptr);
}

TEST_CASE("feature indexes are deterministic across thread counts") {
  SourceDateEpoch epoch;
  tests::TempDir temp;
  std::map<std::string, std::string> first;
  for (unsigned threads : {1u, 2u, 8u}) {
    BuildOptions options;
    options.layout = Layout::Both;
    options.groups = {covidx::Group{"trio", {"S1", "S2"}}};
    options.threads = threads;
    CovIndex index;
    index.build(tiny_bam().string(), options);
    const fs::path out = temp.path() / std::to_string(threads);
    index.save(out.string(), "tiny");
    const auto files = snapshot(out / "tiny.covidx");
    if (first.empty()) {
      first = files;
    } else {
      CHECK(files == first);
    }
  }
}

TEST_CASE("feature query validation") {
  CovIndex empty;
  CHECK_THROWS_AS((void)empty.feature_count("S1", "chr1", 0, 10),
                  covidx::QueryError);

  const CovIndex index = build_tiny();
  CHECK(index.build_stats().features == 12);
  CHECK_THROWS_AS((void)index.feature_count("nope", "chr1", 0, 10),
                  covidx::QueryError);
  CHECK_THROWS_AS((void)index.feature_count("S1", "nope", 0, 10),
                  covidx::QueryError);
  CHECK_THROWS_AS((void)index.feature_count("S1", "chr1", 5, 5),
                  covidx::QueryError);
  CHECK_THROWS_AS((void)index.feature_count("S1", "chr1", 0, 1001),
                  covidx::QueryError);
}
