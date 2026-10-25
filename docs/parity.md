# covidx parity

This document is normative for parity testing.

## Definition

For a window `[X,Y)`, `covidx.covered_count` is the number of reference positions covered by at least one accepted alignment span under the configured policy. That is exactly default `bedtools coverage` output column 5, `num_bases_covered`:

```
chrom  start  end  num_features  num_bases_covered  length  fraction
```

It is **not** `bedtools coverage -counts` (column 4), which counts how many B features overlap the window.

## Covidx policy

| Operation | Reference span | Covered |
|---|---:|---|
| `M`, `=`, `X` | `len` | all; only A/C/G/T under `--ambiguous-read-bases exclude` |
| `D` | `len` | all if `--count-deletions` (default), otherwise a gap |
| `N` | `len` | never; splits blocks |
| `I`, `S`, `H`, `P` | 0 | ignored; do not split |
| `B` | undefined | `InputError` |

Default flag filtering excludes `0x704` (unmapped, secondary, QC-fail, duplicate). Supplementary `0x800` is excluded by default independently of `exclude_flags`. MAPQ 255 is treated as “unavailable” and passes any threshold.

## `samtools depth` mapping

Let:

```
D       = 0x704                         # samtools default filter-out
E       = covidx effective exclude mask
        = (exclude_flags | 0x800) & ~(include_supplementary ? 0x800 : 0)
remove  = D & ~E                        # pass to samtools -g
add     = E & ~D                        # pass to samtools -G
```

Run:

```bash
samtools depth -aa -q 0 -Q <min_mapq> \
  -g <remove> -G <add> \
  [-J if count_deletions] \
  -r <contig> in.bam
```

Notes:

- `samtools depth -G` adds to the default filter-out list; `-g` removes flags from that list. `--excl-flags 0` alone does not remove the defaults.
- `-J` includes deletions; omit it for `--no-count-deletions`.
- `-aa` outputs positions with zero depth and unused reference sequences.
- `-q 0` matches covidx's lack of a base-quality filter.
- `-s` overlap removal changes depth values, not covered/not-covered, so it is optional for coverage parity.
- `samtools depth` includes supplementary alignments by default; covidx does not. Use `-G 0x800` for the default covidx policy.
- `samtools depth` is not a valid reference for `--ambiguous-read-bases exclude`; use the naive checker or a custom depth script that inspects `bam_seqi`.

Policy matrix (derive flags from the formula, do not copy blindly):

| Covidx policy | Effective `E` | Samtools flags |
|---|---:|---|
| defaults | `0xF04` | `-G 0x800 -J` |
| `--include-supplementary` | `0x704` | `-J` |
| `--no-count-deletions` | `0xF04` | `-G 0x800` |
| `--min-mapq 20` | `0xF04` | `-G 0x800 -Q 20 -J` |
| `--exclude-flags 0` | `0x800` | `-g 0x704 -G 0x800 -J` |
| `--exclude-flags 0x604` | `0xE04` | `-g 0x100 -G 0x800 -J` |
| `--ambiguous-read-bases cover` | any | valid |
| `--ambiguous-read-bases exclude` | any | not valid; use `data-raw/naive_reference.py` |

## `bedtools coverage` procedure

1. Build `blocks.bed` from the same BAM with the same policy: accepted alignments only; `M/=/X` contribute; `D` contributes or splits; `N` splits; `I/S/H/P` do not split; one 0-based half-open interval per contributing block.
2. Build `windows.bed` with the same windows used for `covidx-query` (BED3).
3. Run:

   ```bash
   bedtools coverage -a windows.bed -b blocks.bed > coverage.tsv
   ```

4. Compare `covidx.covered_count` to column 5 for every row. Also compare `uncovered = length - covered` and `covered_fraction` to column 7 at six decimals.

Do not compare to `-counts` or column 4. Do not pass the BAM directly as B unless BAM/`-split` interpretation has been verified against the CIGAR policy; the normative comparison uses `blocks.bed`.

`data-raw/check_bedtools_parity.py` implements this procedure.

## Scope

For `--regions` builds, generate blocks only inside the scope and compare windows inside the same scope. Queries that leave a scope interval are `QueryError`; do not silently compare them as zero coverage.

## Tools

- `samtools` 1.21 is the canonical CI version; 1.16+ is supported.
- `bedtools` 2.31.1 is the canonical development version; 2.30+ is supported.
- `data-raw/check_samtools_parity.py`, `data-raw/check_bedtools_parity.py`, and `data-raw/naive_reference.py` are the developer harnesses.
- If a tool is missing, parity tests must skip with a clear message, not pass silently.
