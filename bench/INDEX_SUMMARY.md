# Format 1.1.0 index summary

Measured 2026-10-25 on Debian trixie, 32 vCPU, 117 GiB RAM, GCC 14.2.0, htslib 1.21, Release build. Indexes were built with `SOURCE_DATE_EPOCH=0`.

`covidx` format 1.1.0 stores two payloads in every unit file: canonical coverage runs and one feature span per accepted alignment. The coverage payload is unchanged from 1.0.0. The feature payload is new. `bench/index_summary.py` reads the manifest and unit headers and splits each unit file into coverage bytes, feature bytes, and fixed overhead, without reading the BAM.

Reproduce either table with:

```bash
SOURCE_DATE_EPOCH=0 ./build-release/covidx-build \
  --bam data-raw/multi-sample.bam --out /tmp/covidx-summary \
  --layout both --threads 8
python3 bench/index_summary.py /tmp/covidx-summary/multi-sample.covidx
```

## Five-sample BAM, `--layout both`

Source: `data-raw/multi-sample.bam`, 847.4 MiB, 20,456,242 accepted alignments, 263 contigs.

```text
format 1.1.0  builder 0.2.0  layout both
index: 263 contigs  6 units  7,075,359 runs  40,912,484 logical features (20,456,242 stored)
payload: 187.97 MiB (188.03 MiB on disk)
breakdown: coverage 127.17 MiB, features 60.73 MiB, overhead 0.07 MiB
```

| Unit | Kind | File MiB | Runs | Logical features | Stored features | Coverage MiB | Feature MiB | B/run | B/feature |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Sq54 | sample | 26.25 | 875,555 | 3,534,590 | 3,534,590 | 15.55 | 10.69 | 18.62 | 3.17 |
| Sq81 | sample | 18.50 | 644,787 | 2,244,441 | 2,244,441 | 11.44 | 7.05 | 18.61 | 3.29 |
| Sq84 | sample | 23.90 | 811,579 | 3,082,440 | 3,082,440 | 14.44 | 9.45 | 18.65 | 3.21 |
| Sq85 | sample | 41.52 | 1,258,067 | 6,667,358 | 6,667,358 | 22.51 | 19.00 | 18.76 | 2.99 |
| Sq99 | sample | 34.96 | 1,123,835 | 4,927,413 | 4,927,413 | 20.40 | 14.55 | 19.03 | 3.10 |
| pooled | pooled | 42.85 | 2,361,536 | 20,456,242 | 0 | 42.83 | 0.00 | 19.02 | 0.00 |

The pooled unit is derived from the five sample units in a `both` build. Its manifest `total_features` is the sum of member sample counts, so summing `units[].total_features` gives 40,912,484 logical features even though the unit files store 20,456,242. The summary tool reports both totals so the distinction is explicit.

## Single-sample BAM, `--layout pooled`

Source: `data-raw/cs_t2t_se_cereale.bam`, 174.8 MiB, 3,534,590 accepted alignments, 263 contigs.

```text
format 1.1.0  builder 0.2.0  layout pooled
index: 263 contigs  1 units  875,555 runs  3,534,590 logical features (3,534,590 stored)
payload: 26.25 MiB (26.30 MiB on disk)
breakdown: coverage 15.55 MiB, features 10.69 MiB, overhead 0.01 MiB
```

A pooled-only build stores feature segments in the pooled unit itself, so logical and stored counts agree.

## What the numbers mean

- The coverage portion of the five-sample index is 127.17 MiB, matching the 127.3 MiB measured for the same BAM under format 1.0.0. Feature storage accounts for the entire increase to 188.03 MiB, about 47.8%.
- Coverage costs about 18.85 bytes per run. This is close to the 20-byte model of two `uint32` positions plus the amortized prefix and bucket entries.
- Stored features cost about 3.11 bytes each on 150 bp reads. Each feature contributes a sorted start and a sorted end; most deltas fit in one or two LEB128 bytes, plus 12 bytes of block indexing per 64 values.
- Feature size follows the number of accepted alignments and the coordinate distribution, not covered breadth alone. A saturated region can have very few coverage runs and millions of features. The stress matrix shows four million 150 bp reads producing six runs but a 9.4 MiB index.
- Group and `both`-layout pooled units do not duplicate feature segments. They sum member sample counts at query time. A pooled-only build stores the features once in the pooled unit.
- Headers, 48-byte contig metadata records, alignment padding, TSVs, and the manifest are small for normal indexes: 77,212 bytes for the six-unit, 263-contig example. Per-contig overhead dominates only for indexes with thousands of tiny contigs.
- The build/query optimization pass did not change the format or index size. See `bench/STRESS_RESULTS.md` for the before/after measurements.
