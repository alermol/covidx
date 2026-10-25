#include "doctest.h"

#include <cstdint>
#include <filesystem>
#include <random>

#include "covidx/cov_index.h"
#include "test_util.h"

namespace fs = std::filesystem;
using covidx::CovIndex;
using covidx::OpenOptions;

TEST_CASE("mmap and heap loads return identical answers") {
  tests::TempDir temp;
  {
    CovIndex index;
    index.build((fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam").string());
    index.save(temp.path().string(), "tiny");
  }
  const fs::path directory = temp.path() / "tiny.covidx";
  OpenOptions mmap_options;
  mmap_options.mmap = true;
  OpenOptions heap_options;
  heap_options.mmap = false;

  CovIndex mmap_index;
  mmap_index.load(directory.string(), mmap_options);
  CovIndex heap_index;
  heap_index.load(directory.string(), heap_options);

  std::mt19937 rng(123);
  for (int i = 0; i < 2000; ++i) {
    const std::string unit = i % 2 == 0 ? "S1" : "S2";
    const uint64_t x = rng() % 1000;
    const uint64_t y = x + 1 + rng() % (1000 - x);
    CHECK(mmap_index.covered_count(unit, "chr1", x, y) ==
          heap_index.covered_count(unit, "chr1", x, y));
    CHECK(mmap_index.covered_fraction(unit, "chr1", x, y) ==
          heap_index.covered_fraction(unit, "chr1", x, y));
  }
  CHECK(mmap_index.total_covered("S2", "chr1") ==
        heap_index.total_covered("S2", "chr1"));
}
