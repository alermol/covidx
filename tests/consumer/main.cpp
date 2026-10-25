#include <cstdio>

#include "covidx/cov_index.h"

int main() {
  covidx::CovIndex index;
  const covidx::BuildStats& stats = index.build_stats();
  std::printf("covidx consumer ok: examined=%llu\n",
              static_cast<unsigned long long>(stats.records_examined));
  return 0;
}
