#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>

#include "covidx/cov_index.h"

namespace {

volatile uint64_t sink = 0;

void measure(const covidx::CovIndex& index, const std::string& unit,
             const std::string& contig, uint64_t length, uint64_t window,
             int iterations) {
  std::mt19937_64 rng(12345 + window);
  const uint64_t span = length > window ? length - window + 1 : 1;
  const auto started = std::chrono::steady_clock::now();
  uint64_t sum = 0;
  for (int i = 0; i < iterations; ++i) {
    const uint64_t start = (window >= length) ? 0 : rng() % span;
    sum += index.covered_count(unit, contig, start, start + window);
  }
  const auto elapsed = std::chrono::steady_clock::now() - started;
  const double ns =
      std::chrono::duration<double, std::nano>(elapsed).count() / iterations;
  sink += sum;
  std::printf("%llu\t%d\t%.2f\t%llu\n", static_cast<unsigned long long>(window),
              iterations, ns, static_cast<unsigned long long>(sum));
}

} // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr,
                 "usage: covidx_bench_query <index.covidx> [unit] [contig]\n");
    return 2;
  }
  covidx::CovIndex index;
  index.load(argv[1]);
  const std::string unit = argc > 2 ? argv[2] : index.units().front().name;
  std::string contig;
  uint64_t length = 0;
  for (const covidx::ContigInfo& item : index.contigs()) {
    if (item.length > 0) {
      contig = item.name;
      length = item.length;
      break;
    }
  }
  if (contig.empty()) {
    std::fprintf(stderr, "index has no non-empty contig\n");
    return 2;
  }
  if (argc > 3) {
    contig = argv[3];
    const covidx::ContigInfo* info = index.contig(contig);
    if (info == nullptr) {
      std::fprintf(stderr, "unknown contig %s\n", contig.c_str());
      return 2;
    }
    length = info->length;
  }

  std::printf("window\titerations\tns_per_query\tchecksum\n");
  for (uint64_t window : {uint64_t{1}, uint64_t{100}, uint64_t{10000},
                          uint64_t{1000000}, length}) {
    if (window == 0)
      continue;
    const uint64_t effective = window > length ? length : window;
    measure(index, unit, contig, length, effective, 200000);
  }
  return 0;
}
