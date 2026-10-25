# covidx

[![CI](https://github.com/alermol/covidx/actions/workflows/ci.yml/badge.svg)](https://github.com/alermol/covidx/actions/workflows/ci.yml) [![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

`covidx` builds a compact index from a coordinate-sorted BAM and answers exact breadth-of-coverage queries without rescanning the BAM. It matches default `bedtools coverage` column 5 (`num_bases_covered`) and stores per-alignment spans for column 4 (`num_features`, equivalent to `bedtools coverage -counts`). Use covidx for repeated panel, exome, or QC breadth queries; use bedtools or samtools depth for one-off work or per-base depth.

## Build

Requires a C++23 compiler, CMake 3.20+, htslib 1.16+, zlib, and a POSIX system.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCOVIDX_BUILD_TESTS=OFF
cmake --build build -j
cmake --install build --prefix "$HOME/.local"
```

Add `-G Ninja` if Ninja is installed. If htslib is not in the default pkg-config path, set `PKG_CONFIG_PATH` to its `.pc` directory.

## Quickstart

```bash
./build/covidx-build --bam tests/data/tiny.bam --out indexes
./build/covidx-query --index indexes/tiny.covidx --list-units
./build/covidx-query --index indexes/tiny.covidx --unit S1 --region chr1:1-1000
```

```text
#chrom	start	end	unit	num_features	covered	uncovered	covered_fraction
chr1	0	1000	S1	6	60	940	0.060000
```

`--region` is 1-based inclusive; output coordinates are always 0-based half-open. `num_features` is bedtools coverage column 4 and `covered` is column 5. `covered_fraction = covered / (end - start)`, six decimals.

## Queries

```bash
# BED3+ (0-based half-open); extra columns ignored
covidx-query --index indexes/cohort.covidx --unit S1 --bed targets.bed

# Region (1-based inclusive); add --zero-based for 0-based half-open
covidx-query --index indexes/cohort.covidx --unit S1 --region chr1:1000-2000

# Tile every indexed scope interval
covidx-query --index indexes/cohort.covidx --unit S1 --window 10000 --step 5000

# Inspect the index
covidx-query --index indexes/cohort.covidx --list-units
covidx-query --index indexes/cohort.covidx --list-contigs
```

Every query row includes `num_features` before `covered`; there is no opt-in flag. `--unit` is optional when the index has exactly one unit. Use `--bed -` for stdin, `--output FILE` for a file, and `--no-header` to omit the header. BED input is 0-based half-open, extra columns are ignored, and output preserves input order. `--step` defaults to the window size. Exit codes: `0` success, `1` usage, `2` input/query, `3` output, `4` internal.

## Build an index

```bash
covidx-build --bam cohort.bam --out indexes
covidx-build --bam cohort.bam --out indexes --layout both \
  --group trio=NA12878,NA12891,NA12892
covidx-build --bam cohort.bam --out indexes --layout pooled
```

- BAM only, coordinate-sorted with `.bai`/`.csi`; contig lengths ≤ `INT32_MAX`.
- `--layout per-sample` (default), `pooled`, or `both`.
- The default layout needs an `RG` tag with a declared `SM`. Reads with no or unknown `RG` are an error unless you pass `--unknown-sample NAME` or use `--layout pooled`.
- Add fixed groups with `--group NAME=S1,S2` or `--groups-file FILE` (TSV).
- Restrict scope with `--contigs chr1,chr2`, `--contigs-file FILE`, or `--regions BED` (0-based half-open). `--regions` cannot be combined with the contig options.
- Every index stores one feature per accepted alignment over its full reference span (`M/D/N/=/X`); zero-reference-span alignments have no feature. This provides the mandatory `num_features` query column and format `1.1.0`.
- `--overwrite` atomically replaces an existing index.
- `--threads N` sets workers (default: up to 8). Build parallelism is per contig, so a single-contig BAM will not use more than one worker.

Queries must stay inside one indexed scope interval; out-of-scope windows are an error, not a zero.

### Policy

| Option | Default | Effect |
|---|---|---|
| `--min-mapq N` | 0 | Skip reads with MAPQ < N |
| `--exclude-flags MASK` | `0x704` | Skip matching records |
| `--include-supplementary` | off | Include `0x800` alignments |
| `--count-deletions` | on | `D` covered; off treats it as a gap |
| `--ambiguous-read-bases` | `cover` | `cover` or `exclude` ambiguous bases |

`M`, `=`, and `X` contribute; under `--ambiguous-read-bases exclude`, only A/C/G/T read bases count. `D` is covered by default and a gap with `--no-count-deletions`. `N` always splits blocks and is never covered. `I`, `S`, `H`, and `P` are ignored and do not split.

Default filtering removes unmapped, secondary, QC-fail, and duplicate records; supplementary alignments are excluded separately. MAPQ 255 passes any cutoff. The policy is recorded in `manifest.json`; changing it requires a rebuild. Set `SOURCE_DATE_EPOCH` for byte-identical builds across thread counts.

## Exactness

For the same accepted blocks, `covered` equals default `bedtools coverage` column 5. Raw-BAM comparisons need `-split` plus matching filters; `--no-count-deletions` and `--ambiguous-read-bases exclude` require BED blocks. `num_features` equals default `bedtools coverage` column 4 for the same accepted alignments; the full span is independent of deletion and ambiguous-base coverage policy. See `docs/parity.md`.

## Library

```cpp
#include <cstdint>
#include <covidx/cov_index.h>
#include <iostream>

int main() {
  covidx::CovIndex index;
  index.load("indexes/cohort.covidx");
  std::cout << index.covered_count("S1", "chr1", 1000, 2000) << '\n';
  std::cout << index.feature_count("S1", "chr1", 1000, 2000) << '\n';
}
```

```cmake
find_package(covidx CONFIG REQUIRED)
target_link_libraries(app PRIVATE covidx::covidx)
```

If covidx is installed outside a default prefix, add `-DCMAKE_PREFIX_PATH="$HOME/.local"` when configuring your project.

`CovIndex` is immutable after build/load and safe for concurrent queries. `covered_fraction()`, `units()`, `contigs()`, and `scope()` are also available. For many windows on one contig, `feature_counts()` returns a vector of counts in one feature-array pass; `covidx-query` uses it automatically for large BED and tiling inputs.

## Format and performance

An index is `<out>/<name>.covidx/` with `manifest.json`, TSVs, and per-unit files containing canonical runs, prefix sums, bucket tables, and block-varint alignment spans. Format `1.1.0` is the only format read or written. Indexes from earlier releases must be rebuilt. Loads validate CRCs and structure by default; saves are atomic. See `docs/format.md`.

Index size follows the number of covered runs plus one compressed span per accepted alignment, not covered breadth alone. Queries are independent of BAM size and window length. A five-sample, 847 MiB BAM with 20.5M alignments produced a 188 MiB index: 127 MiB of coverage runs, 61 MiB of feature spans (about 3.1 bytes per stored feature), and 0.07 MiB of headers and metadata. Large BED and tiling queries use a batched feature path; on this index, 500k BED rows take about 1.5 s and a 604k-row 1 kb tiling takes about 1.2 s. Run `python3 bench/index_summary.py <index>.covidx` to break down any index. See `bench/INDEX_SUMMARY.md` and `bench/STRESS_RESULTS.md`.

Limits: breadth only (no per-base depth), BAM only, fixed build policy.

## Development

```bash
cmake -S . -B build -DCOVIDX_BUILD_TESTS=ON -DCOVIDX_WERROR=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Parity tests use samtools and bedtools when available. See `CHANGELOG.md`, `THIRD_PARTY_NOTICES.md`, and `data-raw/`.

## License

MIT. See [`LICENSE`](LICENSE).
