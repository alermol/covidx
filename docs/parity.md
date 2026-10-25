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

`data-raw/check_bedtools_parity.py` implements this procedure. It validates column 5 only. Column 4 is not part of format 1.0.0: it requires per-alignment spans, not the merged run table.

## Real-data validation

The 20 kb `CBFFKZ...` contigs used in earlier checks, such as `CBFFKZ010000001.1`, are unplaced sequences with zero alignments. They are useful for zero-handling tests, but a parity run on such a contig is vacuous and must not be reported as real-data evidence. The harnesses reject a run whose selected windows contain no covered positions unless `--allow-empty` is explicit.

Use a read-rich region on a large contig. For `data-raw/cs_t2t_se_cereale.bam`, `CP155612.1:133001-143000` is a 10 kb region with 1,102 covered positions and 42 alignments. A scoped build for this region is small:

```bash
printf 'CP155612.1\t133000\t143000\n' > region.bed
SOURCE_DATE_EPOCH=0 covidx-build \
  --bam data-raw/cs_t2t_se_cereale.bam \
  --out /tmp/covidx-real --regions region.bed --name real
```

Then run either harness with the same region:

```bash
python3 data-raw/check_samtools_parity.py \
  --bam data-raw/cs_t2t_se_cereale.bam \
  --contig CP155612.1 --region CP155612.1:133001-143000 \
  --index /tmp/covidx-real/real.covidx --unit Sq54

python3 data-raw/check_bedtools_parity.py \
  --bam data-raw/cs_t2t_se_cereale.bam \
  --contig CP155612.1 --region CP155612.1:133001-143000 \
  --index /tmp/covidx-real/real.covidx --unit Sq54
```

`--region` is 1-based inclusive, matching `covidx-query --region`. It generates windows of 1..50 bp inside the region and restricts the reference scan to a BED of those windows. For a larger or irregular window set, pass `--windows-bed` with BED3 windows (0-based half-open). The samtools and bedtools harnesses use `samtools depth -b` without `-aa` in bounded mode, so they do not scan or allocate the whole contig.

For a per-sample unit in a multi-sample BAM, pass `--sample <name>`. The harness reads `<index>/samples.tsv` and filters the reference BAM to that sample's read groups. `--read-group <id[,id...]>` is the explicit form and handles group units. Without this filter, `samtools depth` includes every sample and per-sample parity is invalid.

Both supplied BAMs contain only primary, non-duplicate, non-supplementary reads with MAPQ >= 20. The supplementary, secondary, duplicate, QC-fail, and `--min-mapq 20` policy rows cannot be distinguished on these files. Those rows are covered by synthetic unit tests; real-data runs must state which policies were actually exercised.

The harnesses check the index manifest policy against the requested policy before comparing. A mismatched index is an error, not a silent comparison.

For `--ambiguous-read-bases exclude`, the index and the naive checker must use the same policy. `tests/data/tiny.bam` contains `r_star` with `SEQ=*`, which the exclude policy rejects. Filtering out that record leaves `r_ambig`, whose sequence contains two `N` bases; on `chr1` for sample `S1`, the region `chr1:1-1000` then has 60 covered positions under `cover` and 58 under `exclude`. This exercises the ambiguous-base path without a real reference genome.

## Scope

For `--regions` builds, generate blocks only inside the scope and compare windows inside the same scope. Queries that leave a scope interval are `QueryError`; do not silently compare them as zero coverage.

## Tools

- `samtools` 1.21 is the canonical CI version; 1.16+ is supported.
- `bedtools` 2.31.1 is the canonical development version; 2.30+ is supported.
- `data-raw/check_samtools_parity.py`, `data-raw/check_bedtools_parity.py`, and `data-raw/naive_reference.py` are the developer harnesses.
- The harnesses accept `--region`, `--windows-bed`, `--sample`, `--read-group`, and `--allow-empty`; each verifies the index policy in `manifest.json` before comparing.
- `check_samtools_parity.py` and `check_bedtools_parity.py` reject `--ambiguous-read-bases exclude`; use `naive_reference.py` for that policy.
- If a tool is missing, parity tests must skip with a clear message, not pass silently.
