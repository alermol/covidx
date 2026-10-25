#include "doctest.h"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <string>
#include <thread>
#include <vector>

#include "covidx/cov_index.h"
#include "test_util.h"

namespace fs = std::filesystem;
using covidx::CovIndex;

TEST_CASE("concurrent read-only queries return identical answers") {
  tests::TempDir temp;
  {
    CovIndex index;
    index.build((fs::path(COVIDX_TEST_DATA_DIR) / "tiny.bam").string());
    index.save(temp.path().string(), "tiny");
  }
  CovIndex index;
  index.load((temp.path() / "tiny.covidx").string());

  std::vector<std::pair<std::string, uint64_t>> expected;
  for (const char* unit : {"S1", "S2", "S3", "S4"}) {
    for (uint64_t x = 0; x < 1000; x += 37) {
      expected.emplace_back(
          unit, index.covered_count(unit, "chr1", x,
                                    std::min<uint64_t>(1000, x + 83)));
    }
  }

  std::atomic<bool> failed{false};
  std::vector<std::thread> threads;
  for (int thread = 0; thread < 8; ++thread) {
    threads.emplace_back([&] {
      try {
        for (int round = 0; round < 200; ++round) {
          std::size_t i = 0;
          for (const char* unit : {"S1", "S2", "S3", "S4"}) {
            for (uint64_t x = 0; x < 1000; x += 37) {
              const uint64_t got = index.covered_count(
                  unit, "chr1", x, std::min<uint64_t>(1000, x + 83));
              if (got != expected[i].second) {
                failed.store(true);
                return;
              }
              ++i;
            }
          }
        }
      } catch (...) {
        failed.store(true);
      }
    });
  }
  for (auto& thread : threads)
    thread.join();
  CHECK_FALSE(failed.load());
}
