#include "doctest.h"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include <htslib/sam.h>

#include "covidx/detail/bam_reader.h"
#include "covidx/error.h"
#include "test_util.h"

namespace fs = std::filesystem;
using covidx::detail::BamHeaderInfo;
using covidx::detail::BamRegionReader;
using covidx::detail::read_bam_header;

namespace {

fs::path write_bam(const fs::path& dir, std::string_view name,
                   std::string_view header,
                   std::vector<tests::ReadSpec> reads = {}) {
  const fs::path path = dir / name;
  tests::make_bam(path, header, std::move(reads));
  return path;
}

tests::ReadSpec one_read(std::string rg = "rg1") {
  return tests::ReadSpec{
      "read1",      0, 0, 10, 60, {{'M', 10}}, tests::repeat_sequence(10),
      std::move(rg)};
}

} // namespace

TEST_CASE("reads a valid tiny BAM header and sample map") {
  const fs::path path = fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam";
  const BamHeaderInfo info = read_bam_header(path.string());
  REQUIRE(info.contigs.size() == 2);
  CHECK(info.contigs[0].name == "chr1");
  CHECK(info.contigs[0].length == 1000);
  CHECK(info.contigs[1].name == "chr2");
  CHECK(info.contigs[1].length == 500);
  REQUIRE(info.samples.size() == 4);
  CHECK(info.samples[0].name == "S1");
  CHECK(info.samples[1].name == "S2");
  CHECK(info.samples[2].name == "S3");
  CHECK(info.samples[3].name == "S4");
  REQUIRE(info.samples[1].read_groups.size() == 2);
  CHECK(info.samples[1].read_groups[0] == "rg2");
  CHECK(info.samples[1].read_groups[1] == "rg3");
  CHECK(info.samples[2].read_groups.size() == 1);
  CHECK(info.rg_to_sample.at("rg2") == 1);
  CHECK(info.rg_to_sample.at("rg3") == 1);
  CHECK(info.rg_to_sample.at("rg4") == 2);
  CHECK(info.header_so == "coordinate");
  CHECK(info.file_size > 0);
  CHECK(fs::exists(info.index_path));
}

TEST_CASE("header validation errors are distinct") {
  tests::TempDir temp;
  CHECK_THROWS_AS(
      (void)read_bam_header(
          write_bam(temp.path(), "no_hd.bam",
                    "@SQ\tSN:chr1\tLN:100\n@RG\tID:rg1\tSM:S1\n", {one_read()})
              .string()),
      covidx::InputError);
  CHECK_THROWS_AS(
      (void)read_bam_header(
          write_bam(temp.path(), "queryname.bam",
                    "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:100\n"
                    "@RG\tID:rg1\tSM:S1\n",
                    {one_read()})
              .string()),
      covidx::InputError);
  CHECK_THROWS_AS(
      (void)read_bam_header(
          write_bam(temp.path(), "zero_ln.bam",
                    "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:0\n", {})
              .string()),
      covidx::InputError);
  CHECK_THROWS_AS(
      (void)read_bam_header(
          write_bam(temp.path(), "long_ln.bam",
                    "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:3000000000\n",
                    {})
              .string()),
      covidx::InputError);
  CHECK_THROWS_AS(
      (void)read_bam_header(
          write_bam(temp.path(), "dup_rg.bam",
                    "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n"
                    "@RG\tID:rg1\tSM:S1\n@RG\tID:rg1\tSM:S2\n",
                    {one_read()})
              .string()),
      covidx::InputError);
  CHECK_THROWS_AS(
      (void)read_bam_header(
          write_bam(temp.path(), "no_sm.bam",
                    "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n"
                    "@RG\tID:rg1\n",
                    {one_read()})
              .string()),
      covidx::InputError);
}

TEST_CASE("missing index and non-BAM inputs are rejected") {
  tests::TempDir temp;
  const fs::path bam =
      write_bam(temp.path(), "no_index.bam", tests::kTinyHeader, {one_read()});
  fs::remove(bam.string() + ".bai");
  CHECK_THROWS_AS((void)read_bam_header(bam.string(), true),
                  covidx::InputError);
  CHECK_NOTHROW((void)read_bam_header(bam.string(), false));

  const fs::path sam = temp.path() / "input.sam";
  tests::write_file(sam, "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n");
  CHECK_THROWS_AS((void)read_bam_header(sam.string(), false),
                  covidx::InputError);
  CHECK_THROWS_AS((void)read_bam_header("https://example.com/x.bam", false),
                  covidx::InputError);
  CHECK_THROWS_AS(
      (void)read_bam_header((temp.path() / "missing.bam").string(), false),
      covidx::InputError);
}

TEST_CASE("sample_id_for_read resolves and rejects RG tags") {
  tests::TempDir temp;
  const fs::path bam = write_bam(temp.path(), "sample.bam", tests::kTinyHeader,
                                 {one_read("rg2")});
  BamHeaderInfo info = read_bam_header(bam.string());

  samFile* fp = sam_open(bam.string().c_str(), "r");
  REQUIRE(fp != nullptr);
  sam_hdr_t* hdr = sam_hdr_read(fp);
  REQUIRE(hdr != nullptr);
  bam1_t* record = bam_init1();
  REQUIRE(sam_read1(fp, hdr, record) >= 0);
  CHECK(covidx::detail::sample_id_for_read(*record, info, "") ==
        info.rg_to_sample.at("rg2"));
  bam_destroy1(record);
  sam_hdr_destroy(hdr);
  sam_close(fp);

  const fs::path unknown = write_bam(temp.path(), "unknown.bam",
                                     tests::kTinyHeader, {one_read("rgX")});
  info = read_bam_header(unknown.string());
  fp = sam_open(unknown.string().c_str(), "r");
  hdr = sam_hdr_read(fp);
  record = bam_init1();
  REQUIRE(sam_read1(fp, hdr, record) >= 0);
  CHECK_THROWS_AS((void)covidx::detail::sample_id_for_read(*record, info, ""),
                  covidx::InputError);
  const uint32_t synthetic = static_cast<uint32_t>(info.samples.size());
  info.samples.push_back(covidx::detail::SampleRecord{"UNK", {}});
  info.sample_by_name.emplace("UNK", synthetic);
  CHECK(covidx::detail::sample_id_for_read(*record, info, "UNK") == synthetic);
  CHECK(covidx::detail::sample_id_for_read(*record, info, "S1") == 0);
  bam_destroy1(record);
  sam_hdr_destroy(hdr);
  sam_close(fp);
}

TEST_CASE("region reader returns coordinate-sorted reads once") {
  const fs::path path = fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam";
  BamRegionReader reader(path.string(), 0, {{0, 1000}});
  bam1_t* record = bam_init1();
  int64_t last = -1;
  int count = 0;
  while (reader.next(*record)) {
    CHECK(record->core.pos >= last);
    last = record->core.pos;
    ++count;
  }
  CHECK(count == static_cast<int>(tests::tiny_reads().size()));
  bam_destroy1(record);

  BamRegionReader empty(path.string(), 1, {{0, 500}});
  bam1_t* other = bam_init1();
  CHECK_FALSE(empty.next(*other));
  bam_destroy1(other);
}

TEST_CASE("region reader deduplicates reads spanning a scope gap") {
  const fs::path path = fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam";
  // r_n spans 300..320 and overlaps both intervals.
  BamRegionReader reader(path.string(), 0, {{300, 310}, {315, 325}});
  bam1_t* record = bam_init1();
  int count = 0;
  while (reader.next(*record))
    ++count;
  CHECK(count == 1);
  bam_destroy1(record);
}

TEST_CASE("region reader restarts coordinate order at each interval") {
  tests::TempDir temp;
  std::vector<tests::ReadSpec> reads = {
      {"span", 0, 0, 100, 60, {{'M', 100}}, tests::repeat_sequence(100), "rg1"},
      {"later", 0, 0, 150, 60, {{'M', 20}}, tests::repeat_sequence(20), "rg1"},
  };
  const fs::path bam =
      write_bam(temp.path(), "spans.bam", tests::kTinyHeader, std::move(reads));
  // The second interval re-emits "span" (POS 100) after the first interval
  // already returned "later" (POS 150). That is valid input: the new iterator
  // restarts at the interval's overlapping prefix.
  BamRegionReader reader(bam.string(), 0, {{160, 161}, {195, 196}});
  bam1_t* record = bam_init1();
  std::vector<int64_t> positions;
  while (reader.next(*record))
    positions.push_back(record->core.pos);
  bam_destroy1(record);
  CHECK(positions == std::vector<int64_t>{100, 150});
}

TEST_CASE("truncated BAM produces an input error") {
  tests::TempDir temp;
  const fs::path source = fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam";
  const fs::path bam = temp.path() / "truncated.bam";
  fs::copy_file(source, bam);
  fs::copy_file(source.string() + ".bai", bam.string() + ".bai");
  fs::resize_file(bam, 350);
  BamRegionReader reader(bam.string(), 0, {{0, 1000}});
  bam1_t* record = bam_init1();
  CHECK_THROWS_AS((void)reader.next(*record), covidx::InputError);
  bam_destroy1(record);
}
