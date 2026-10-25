# Benchmark baseline

Environment: Debian trixie, 32 vCPU, GCC 14.2.0, htslib 1.21, Release build.

Command:

```bash
cmake -S . -B build-release -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCOVIDX_BUILD_BENCH=ON -DCOVIDX_BUILD_TESTS=OFF
cmake --build build-release -j
./build-release/covidx-build --bam tests/data/tiny.bam --out /tmp/covidx-bench
./build-release/covidx_bench_query /tmp/covidx-bench/tiny.covidx S1 chr1
```

The tiny fixture is a constant-time smoke test. Query cost does not grow with window length beyond constant overhead.

The scale build benchmark uses 10 contigs of 2,000,000,000 bp each (20 Gb total) with 1,000,000 runs per contig (10,000,000 runs total):

```bash
./build/covidx_bench_build --contigs 10 --runs-per-contig 1000000 --length 2000000000
# contigs=10 runs=10000000 covered=1000000000 elapsed_s=0.331
# estimated_index_mb=190.7 checksum=497857
```

The index-size model is ~20 bytes per run plus bucket entries; the estimated 190.7 MB is the run portion (10M * 20 B). Format 1.1.0 also stores feature-count segments, so this model covers the coverage-run portion only. For measured feature overhead and a full index breakdown, see `STRESS_RESULTS.md` and run `python3 bench/index_summary.py <index.covidx>`.

| Window | Iterations | ns/query | Notes |
|---:|---:|---:|---|
| 1 | 200000 | 340 | fixture |
| 100 | 200000 | 350 | fixture |
| 10000 | 200000 | 252 | clamped to 1000 |
| 1000000 | 200000 | 252 | clamped to 1000 |
| 1000 | 200000 | 251 | whole contig |

Replace this table with measured values from the release hardware; do not treat the placeholder numbers as a CI gate.
