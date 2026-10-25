#include "doctest.h"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "covidx/detail/bucket_table.h"
#include "covidx/detail/crc32.h"
#include "covidx/detail/unit_file.h"
#include "covidx/error.h"
#include "covidx/format.h"
#include "test_util.h"

namespace fs = std::filesystem;
using covidx::detail::bucket_entries;
using covidx::detail::build_first;
using covidx::detail::choose_w;
using covidx::detail::Run;
using covidx::detail::RunIndexView;
using covidx::detail::UnitContigMeta;
using covidx::detail::UnitFileReader;
using covidx::detail::UnitFileWriter;

namespace {

struct Arrays {
  std::vector<uint32_t> start;
  std::vector<uint32_t> end;
  std::vector<uint64_t> cum;
  std::vector<uint32_t> first;
  uint32_t w = 20;
  uint64_t total = 0;
  uint64_t max_bucket_population = 1;
};

Arrays make_arrays(const std::vector<Run>& runs, uint64_t length) {
  Arrays arrays;
  arrays.w = choose_w(length, runs.size());
  arrays.cum.assign(runs.size() + 1, 0);
  for (std::size_t i = 0; i < runs.size(); ++i) {
    arrays.start.push_back(static_cast<uint32_t>(runs[i].start));
    arrays.end.push_back(static_cast<uint32_t>(runs[i].end));
    arrays.cum[i + 1] = arrays.cum[i] + (runs[i].end - runs[i].start);
  }
  arrays.total = arrays.cum.back();
  arrays.first = build_first(runs, length, arrays.w);
  return arrays;
}

void append(UnitFileWriter& writer, const Arrays& arrays, uint64_t total,
            uint64_t max_pop = 1) {
  UnitContigMeta meta;
  meta.runs = arrays.start.size();
  meta.total_covered = total;
  meta.max_bucket_population = max_pop;
  meta.w = arrays.w;
  writer.append_contig(meta, arrays.start, arrays.end, arrays.cum,
                       arrays.first);
}

void check_view(const RunIndexView& view, const Arrays& arrays,
                uint64_t length) {
  CHECK(view.length == length);
  CHECK(view.runs == arrays.start.size());
  CHECK(view.total_covered == arrays.total);
  CHECK(view.w == arrays.w);
  REQUIRE(view.start.size() == arrays.start.size());
  for (std::size_t i = 0; i < arrays.start.size(); ++i) {
    CHECK(view.start[i] == arrays.start[i]);
    CHECK(view.end[i] == arrays.end[i]);
  }
  if (arrays.start.empty()) {
    CHECK(view.cum.empty());
    CHECK(view.first.empty());
    return;
  }
  REQUIRE(view.cum.size() == arrays.cum.size());
  for (std::size_t i = 0; i < arrays.cum.size(); ++i)
    CHECK(view.cum[i] == arrays.cum[i]);
  REQUIRE(view.first.size() == arrays.first.size());
  for (std::size_t i = 0; i < arrays.first.size(); ++i)
    CHECK(view.first[i] == arrays.first[i]);
}

std::string patch_byte(std::string bytes, std::size_t offset, uint8_t value) {
  REQUIRE(offset < bytes.size());
  bytes[offset] = static_cast<char>(value);
  return bytes;
}

std::string patch_u32(std::string bytes, std::size_t offset, uint32_t value) {
  REQUIRE(offset + 4 <= bytes.size());
  for (std::size_t i = 0; i < 4; ++i)
    bytes[offset + i] = static_cast<char>(value >> (8 * i));
  return bytes;
}

std::string patch_u64(std::string bytes, std::size_t offset, uint64_t value) {
  REQUIRE(offset + 8 <= bytes.size());
  for (std::size_t i = 0; i < 8; ++i)
    bytes[offset + i] = static_cast<char>(value >> (8 * i));
  return bytes;
}

} // namespace

TEST_CASE("unit metadata encodes and decodes exactly 48 bytes") {
  UnitContigMeta meta;
  meta.segment_offset = 0x1122334455667788ULL;
  meta.runs = 42;
  meta.total_covered = 1000;
  meta.max_bucket_population = 7;
  meta.w = 12;
  meta.flags = 3;
  meta.feature_offset = 0x8877665544332211ULL;
  const auto bytes = meta.encode();
  CHECK(bytes.size() == 48);
  const UnitContigMeta decoded = UnitContigMeta::decode(
      std::span<const uint8_t>(bytes.data(), bytes.size()));
  CHECK(decoded.segment_offset == meta.segment_offset);
  CHECK(decoded.runs == meta.runs);
  CHECK(decoded.total_covered == meta.total_covered);
  CHECK(decoded.max_bucket_population == meta.max_bucket_population);
  CHECK(decoded.w == meta.w);
  CHECK(decoded.flags == meta.flags);
  CHECK(decoded.feature_offset == meta.feature_offset);
  for (std::size_t i = 37; i < 40; ++i)
    CHECK(bytes[i] == 0);
}

TEST_CASE("unit file round-trips with mmap and heap readers") {
  tests::TempDir temp;
  const fs::path path = temp.path() / "u000007.bin";
  const uint64_t length0 = 1000;
  const uint64_t length1 = 500;
  const Arrays arrays0 = make_arrays({{0, 5}, {10, 20}, {100, 200}}, length0);
  const Arrays empty = make_arrays({}, length1);
  {
    UnitFileWriter writer(path.string(), 7, 2);
    append(writer, arrays0, arrays0.total, 3);
    append(writer, empty, 0, 0);
    const uint32_t crc = writer.close_and_crc();
    CHECK(crc == covidx::detail::crc32_file(path.string()));
  }

  CHECK(fs::file_size(path) % 8 == 0);
  const std::string sidecar = tests::read_file(path.string() + ".crc32");
  REQUIRE(sidecar.size() == 9);
  CHECK(sidecar[8] == '\n');
  CHECK(sidecar.substr(0, 8) ==
        covidx::detail::crc32_hex(covidx::detail::crc32_file(path.string())));

  for (bool mmap : {true, false}) {
    UnitFileReader reader(path.string(), {length0, length1}, mmap);
    CHECK(reader.used_mmap() == mmap);
    CHECK(reader.header().unit_id == 7);
    CHECK(reader.header().contig_count == 2);
    CHECK(reader.header().meta_record_size == 48);
    CHECK(reader.header().data_offset == 160);
    CHECK(reader.meta(0).segment_offset == 160);
    CHECK(reader.meta(1).segment_offset == 0);
    CHECK(reader.meta(1).w == 20);
    check_view(reader.view(0), arrays0, length0);
    check_view(reader.view(1), empty, length1);
    reader.validate(covidx::Validation::Full);
    CHECK(reader.file_crc32() == covidx::detail::crc32_file(path.string()));
  }
}

TEST_CASE("unit file handles one hundred contigs and dense runs") {
  tests::TempDir temp;
  const fs::path path = temp.path() / "u000000.bin";
  const uint64_t length = 200;
  std::vector<Run> dense;
  for (uint64_t i = 0; i < length; i += 2)
    dense.push_back(Run{i, i + 1});
  const Arrays arrays = make_arrays(dense, length);
  {
    UnitFileWriter writer(path.string(), 0, 100);
    for (int i = 0; i < 100; ++i)
      append(writer, arrays, arrays.total, 2);
    (void)writer.close_and_crc();
  }
  UnitFileReader reader(path.string(), std::vector<uint64_t>(100, length),
                        true);
  for (uint32_t i = 0; i < 100; ++i) {
    CHECK(reader.meta(i).runs == 100);
    check_view(reader.view(i), arrays, length);
  }
  reader.validate(covidx::Validation::Full);
}

TEST_CASE("unit file supports INT32_MAX contig length") {
  tests::TempDir temp;
  const fs::path path = temp.path() / "u000000.bin";
  const uint64_t length = covidx::kMaxContigLength;
  const std::vector<Run> runs = {{length - 1, length}};
  const Arrays arrays = make_arrays(runs, length);
  {
    UnitFileWriter writer(path.string(), 0, 1);
    append(writer, arrays, arrays.total, 1);
    (void)writer.close_and_crc();
  }
  UnitFileReader reader(path.string(), {length}, true);
  const RunIndexView view = reader.view(0);
  CHECK(view.w == 20);
  CHECK(view.first.size() == bucket_entries(length, 20));
  CHECK(covidx::detail::prefix_covered(view, length - 1) == 0);
  CHECK(covidx::detail::prefix_covered(view, length) == 1);
  reader.validate(covidx::Validation::Full);
}

TEST_CASE("unit file stores and queries multi-block feature counts") {
  tests::TempDir temp;
  const fs::path path = temp.path() / "features.bin";
  const uint64_t length = 1000;
  std::vector<uint32_t> starts;
  std::vector<uint32_t> ends;
  for (uint32_t i = 0; i < 150; ++i) {
    starts.push_back(i / 2);
    ends.push_back(i / 2 + 1);
  }
  {
    UnitFileWriter writer(path.string(), 0, 1, true);
    UnitContigMeta meta;
    writer.append_contig(meta, {}, {}, {}, {}, starts, ends);
    (void)writer.close_and_crc();
  }

  for (bool mmap : {true, false}) {
    UnitFileReader reader(path.string(), {length}, mmap);
    CHECK(reader.has_feature_counts());
    CHECK(reader.feature_total(0) == starts.size());
    CHECK(reader.feature_count(0, 0, 1) == 2);
    CHECK(reader.feature_count(0, 0, 2) == 4);
    CHECK(reader.feature_count(0, 1, 2) == 2);
    CHECK(reader.feature_count(0, 0, 1000) == starts.size());
    CHECK(reader.feature_count(0, 40, 45) == 10);
    reader.validate(covidx::Validation::Full);
  }
}

TEST_CASE("unit file interleaves coverage and feature segments") {
  tests::TempDir temp;
  const fs::path path = temp.path() / "mixed.bin";
  const std::vector<Run> runs = {{0, 5}, {10, 20}};
  const Arrays arrays = make_arrays(runs, 1000);
  std::vector<uint32_t> starts;
  std::vector<uint32_t> ends;
  for (uint32_t i = 0; i < 70; ++i) {
    starts.push_back(i * 2);
    ends.push_back(i * 2 + 1);
  }
  {
    UnitFileWriter writer(path.string(), 0, 2, true);
    UnitContigMeta first;
    first.runs = arrays.start.size();
    first.total_covered = arrays.total;
    first.max_bucket_population = 2;
    first.w = arrays.w;
    writer.append_contig(first, arrays.start, arrays.end, arrays.cum,
                         arrays.first, starts, ends);
    UnitContigMeta second;
    const std::array<uint32_t, 2> second_starts = {5, 6};
    const std::array<uint32_t, 2> second_ends = {10, 12};
    writer.append_contig(second, {}, {}, {}, {}, second_starts, second_ends);
    (void)writer.close_and_crc();
  }

  UnitFileReader reader(path.string(), {1000, 500}, true);
  CHECK(reader.has_feature_counts());
  CHECK(reader.feature_total(0) == 70);
  CHECK(reader.feature_total(1) == 2);
  CHECK(reader.feature_count(0, 0, 10) == 5);
  CHECK(reader.feature_count(0, 4, 7) == 2);
  CHECK(reader.feature_count(1, 4, 7) == 2);
  CHECK(reader.view(0).runs == 2);
  CHECK(reader.view(1).runs == 0);
  reader.validate(covidx::Validation::Full);
}

TEST_CASE("unit file rejects feature segment corruption") {
  tests::TempDir temp;
  const fs::path path = temp.path() / "features.bin";
  const uint64_t length = 1000;
  std::vector<uint32_t> starts;
  std::vector<uint32_t> ends;
  for (uint32_t i = 0; i < 150; ++i) {
    starts.push_back(i);
    ends.push_back(i + 2);
  }
  {
    UnitFileWriter writer(path.string(), 0, 1, true);
    UnitContigMeta meta;
    writer.append_contig(meta, {}, {}, {}, {}, starts, ends);
    (void)writer.close_and_crc();
  }
  const std::string good = tests::read_file(path);
  const std::size_t feature_offset = 112;
  REQUIRE(good.size() > feature_offset + 64);

  const auto expect_bad = [&](std::string bytes) {
    const fs::path bad = temp.path() / "bad.bin";
    tests::write_file(bad, bytes);
    CHECK_THROWS_AS(UnitFileReader(bad.string(), {length}, false),
                    covidx::FormatError);
  };

  expect_bad(patch_byte(good, feature_offset, 'X'));
  expect_bad(patch_u32(good, feature_offset + 12, 32));
  expect_bad(patch_u64(good, feature_offset + 16, 0));
  expect_bad(patch_u64(good, feature_offset + 24, 0));
  expect_bad(patch_u64(good, feature_offset + 56, 1));
  expect_bad(good.substr(0, feature_offset + 80));

  // The block "last" summaries are checked while decoding a query block.
  const std::size_t blocks = (starts.size() + 63) / 64;
  const std::size_t starts_last = feature_offset + 64 + blocks * 8;
  const std::string bad_last = patch_u32(good, starts_last, 70);
  const fs::path bad_path = temp.path() / "bad_last.bin";
  tests::write_file(bad_path, bad_last);
  UnitFileReader reader(bad_path.string(), {length}, false,
                        covidx::Validation::None);
  CHECK_THROWS_AS((void)reader.feature_count(0, 64, 65), covidx::FormatError);
}

TEST_CASE("unit file rejects structural corruption") {
  tests::TempDir temp;
  const fs::path path = temp.path() / "u000000.bin";
  const uint64_t length = 1000;
  const Arrays arrays = make_arrays({{0, 5}, {10, 20}}, length);
  {
    UnitFileWriter writer(path.string(), 3, 1);
    append(writer, arrays, arrays.total, 1);
    (void)writer.close_and_crc();
  }
  const std::string good = tests::read_file(path);

  const auto expect_bad = [&](std::string bytes) {
    const fs::path bad = temp.path() / "bad.bin";
    tests::write_file(bad, bytes);
    CHECK_THROWS_AS(UnitFileReader(bad.string(), {length}, false),
                    covidx::FormatError);
  };

  expect_bad(patch_byte(good, 0, 'X'));
  expect_bad(patch_u32(good, 12, 2));
  expect_bad(patch_u32(good, 16, 40));
  expect_bad(patch_u64(good, 64, 65));          // unaligned segment
  expect_bad(patch_u64(good, 64 + 8, 1000000)); // runs out of range
  expect_bad(patch_u64(good, 64 + 16, 1001));   // total beyond length
  expect_bad(patch_u32(good, 64 + 32, 5));      // w too small
  expect_bad(good.substr(0, good.size() - 8));  // truncated
}

TEST_CASE("unit file full validation catches unsampled array corruption") {
  tests::TempDir temp;
  const fs::path path = temp.path() / "u000000.bin";
  const uint64_t length = 1000;
  std::vector<Run> runs;
  for (uint64_t i = 0; i < 100; i += 2)
    runs.push_back(Run{i, i + 1});
  const Arrays arrays = make_arrays(runs, length);
  {
    UnitFileWriter writer(path.string(), 0, 1);
    append(writer, arrays, arrays.total, 2);
    (void)writer.close_and_crc();
  }
  const std::string good = tests::read_file(path);
  const std::size_t segment = 64 + 48;
  const std::size_t run_count = arrays.start.size();
  const std::size_t cum_offset = segment + 2 * run_count * 4;
  const std::size_t first_offset = cum_offset + (run_count + 1) * 8;

  // Corrupt a middle cum entry; structural sampling does not touch it.
  const fs::path cum_bad = temp.path() / "cum_bad.bin";
  tests::write_file(cum_bad, patch_u64(good, cum_offset + 20 * 8, 999));
  UnitFileReader cum_reader(cum_bad.string(), {length}, false);
  CHECK_THROWS_AS(cum_reader.validate(covidx::Validation::Full),
                  covidx::FormatError);

  // Corrupt a middle bucket entry.
  const fs::path bucket_bad = temp.path() / "bucket_bad.bin";
  tests::write_file(bucket_bad, patch_u32(good, first_offset + 4, 7));
  UnitFileReader bucket_reader(bucket_bad.string(), {length}, false);
  CHECK_THROWS_AS(bucket_reader.validate(covidx::Validation::Full),
                  covidx::FormatError);
}

TEST_CASE("CRC sidecar parsing rejects malformed files") {
  tests::TempDir temp;
  const fs::path sidecar = temp.path() / "x.crc32";
  tests::write_file(sidecar, "cbf43926\n");
  CHECK(covidx::detail::read_crc_sidecar(sidecar.string()) == 0xCBF43926u);

  tests::write_file(sidecar, "CBF43926\n");
  CHECK_THROWS_AS((void)covidx::detail::read_crc_sidecar(sidecar.string()),
                  covidx::FormatError);
  tests::write_file(sidecar, "cbf43926");
  CHECK_THROWS_AS((void)covidx::detail::read_crc_sidecar(sidecar.string()),
                  covidx::FormatError);
  tests::write_file(sidecar, "cbf4392\n");
  CHECK_THROWS_AS((void)covidx::detail::read_crc_sidecar(sidecar.string()),
                  covidx::FormatError);
}
