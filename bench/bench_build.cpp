#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "covidx/detail/bucket_table.h"
#include "covidx/detail/run_index.h"
#include "covidx/format.h"

namespace {

uint64_t parse(const char* text) {
  char* end = nullptr;
  const unsigned long long value = std::strtoull(text, &end, 10);
  if (end == text) {
    std::fprintf(stderr, "invalid integer: %s\n", text);
    std::exit(2);
  }
  return static_cast<uint64_t>(value);
}

} // namespace

int main(int argc, char** argv) {
  uint64_t contigs = 10;
  uint64_t runs_per_contig = 1000000;
  uint64_t length = 2000000000;
  for (int i = 1; i + 1 < argc; i += 2) {
    const std::string option = argv[i];
    if (option == "--contigs")
      contigs = parse(argv[i + 1]);
    else if (option == "--runs-per-contig")
      runs_per_contig = parse(argv[i + 1]);
    else if (option == "--length")
      length = parse(argv[i + 1]);
    else {
      std::fprintf(stderr, "unknown option %s\n", argv[i]);
      return 2;
    }
  }
  if (contigs == 0 || runs_per_contig == 0 || length == 0 ||
      length > covidx::kMaxContigLength) {
    std::fprintf(stderr, "invalid scale parameters\n");
    return 2;
  }

  const auto started = std::chrono::steady_clock::now();
  uint64_t total_runs = 0;
  uint64_t total_covered = 0;
  uint64_t checksum = 0;
  std::mt19937_64 rng(2026);
  for (uint64_t contig = 0; contig < contigs; ++contig) {
    const uint64_t stride = length / (runs_per_contig + 1);
    std::vector<covidx::detail::Run> runs;
    runs.reserve(static_cast<std::size_t>(runs_per_contig));
    for (uint64_t i = 0; i < runs_per_contig; ++i) {
      const uint64_t start = (i + 1) * stride;
      const uint64_t end = start + 100;
      if (end > length)
        break;
      runs.push_back(covidx::detail::Run{start, end});
    }
    covidx::detail::RunIndex index;
    index.reset(length);
    index.set_runs(std::move(runs));
    total_runs += index.runs();
    total_covered += index.total_covered();
    for (int q = 0; q < 1000; ++q) {
      const uint64_t x = rng() % (length - 1000);
      checksum += index.covered_count(x, x + 1000);
    }
  }
  const auto elapsed =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - started);
  const double bytes = static_cast<double>(total_runs) * 20.0;
  std::printf("contigs=%llu runs=%llu covered=%llu elapsed_s=%.3f "
              "estimated_index_mb=%.1f "
              "checksum=%llu\n",
              static_cast<unsigned long long>(contigs),
              static_cast<unsigned long long>(total_runs),
              static_cast<unsigned long long>(total_covered), elapsed.count(),
              bytes / (1024.0 * 1024.0),
              static_cast<unsigned long long>(checksum));
  return 0;
}
