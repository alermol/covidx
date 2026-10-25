# covidx

[![CI](https://github.com/alermol/covidx/actions/workflows/ci.yml/badge.svg)](https://github.com/alermol/covidx/actions/workflows/ci.yml) [![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

`covidx` builds a compact index from a coordinate-sorted BAM and answers exact breadth-of-coverage queries without rescanning the BAM. It matches default `bedtools coverage` column 5 (`num_bases_covered`), not `-counts`. Use covidx for repeated panel, exome, or QC breadth queries; use bedtools or samtools depth for one-off work or per-base depth.

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
#chrom	start	end	unit	covered	uncovered	covered_fraction
chr1	0	1000	S1	60	940	0.060000
```

`--region` is 1-based inclusive; output coordinates are always 0-based half-open. `covered_fraction = covered / (end - start)`, six decimals.

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

`--unit` is optional when the index has exactly one unit. Use `--bed -` for stdin, `--output FILE` for a file, and `--no-header` to omit the header. BED input is 0-based half-open, extra columns are ignored, and output preserves input order. `--step` defaults to the window size. Exit codes: `0` success, `1` usage, `2` input/query, `3` output, `4` internal.

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

For the same accepted blocks, `covered` equals default `bedtools coverage` column 5. Raw-BAM comparisons need `-split` plus matching filters; `--no-count-deletions` and `--ambiguous-read-bases exclude` require BED blocks. See `docs/parity.md`.

## Library

```cpp
#include <cstdint>
#include <covidx/cov_index.h>
#include <iostream>

int main() {
  covidx::CovIndex index;
  index.load("indexes/cohort.covidx");
  std::cout << index.covered_count("S1", "chr1", 1000, 2000) << '\n';
}
```

```cmake
find_package(covidx CONFIG REQUIRED)
target_link_libraries(app PRIVATE covidx::covidx)
```

If covidx is installed outside a default prefix, add `-DCMAKE_PREFIX_PATH="$HOME/.local"` when configuring your project.

`CovIndex` is immutable after build/load and safe for concurrent queries. `covered_fraction()`, `units()`, `contigs()`, and `scope()` are also available.

## Format and performance

An index is `<out>/<name>.covidx/` with `manifest.json`, TSVs, and per-unit files containing canonical runs, prefix sums, and bucket tables. Format `1.0.0` is little-endian; loads validate CRCs and structure by default; saves are atomic. See `docs/format.md`.

Index size follows the number of covered runs, not sequencing depth. Queries are independent of BAM size and window length. See `bench/STRESS_RESULTS.md`.

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
