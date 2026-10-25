# covidx stress-test results

Measured 2026-10-25 on Debian trixie, 32 vCPU, 117 GiB RAM, GCC 14.2.0, htslib 1.21, Release build. BAMs, indexes, and temporary work directories were on tmpfs (`/tmp`). All builds used `SOURCE_DATE_EPOCH=0` unless noted. Every index was built with format `1.1.0`, so every measurement includes mandatory feature-count storage and feature queries.

Reproduce with:

```bash
cmake --build build-release -j
./bench/stress_matrix.py --stage all \
  --work-dir /tmp/covidx-stress \
  --windows 5000 --threads 8
python3 bench/index_summary.py /tmp/covidx-stress/indexes/query_multisample/multi-sample.covidx
```

The harness writes `stress-results.json` plus one TSV per stage. The results from this run are copied to `bench/stress-results/`. `bench/index_summary.py` splits an index into coverage, feature, and metadata bytes; see `bench/INDEX_SUMMARY.md` for the measured format summary.

## Matrix

| Dimension | Values exercised |
|---|---|
| Reads | 10k to 20.46M |
| BAM size | 0.24 MiB to 847 MiB |
| Contigs | 1, 4, 64, 256, 4,096, 8,192; real BAMs have 263 |
| Contig length | 100 kb, 1 Mb, 10 Mb, 500 Mb, 2,147,483,647 bp |
| Read length | 50, 150, 1,000, 10,000 bp |
| CIGAR | `M`, `D`, `N`, `I`, `S`, `H`, `=`/`X`, mixed, 600-op |
| MAPQ | constant, varied, unavailable (255) |
| Flags | clean, duplicate/secondary/supplementary/QC-fail mix |
| Samples | 0, 1, 4, 8, 64, 256 |
| Read groups | declared, unknown, absent |
| Spatial layout | uniform, 90% clustered, unmapped tail |
| Index | BAI, CSI, CSI at INT32_MAX contig length |
| Build | threads 1-30, layout per-sample/pooled/both, MAPQ/flag/indel/base policies, scope, checksum |
| Features | 0 to 20,456,242 stored spans; 40,912,484 logical spans in a `both` layout |
| Query | mmap/heap, CRC on/off, validation none/structural/full, BED 1k-500k rows, 1 kb-100 Mb tiling, feature counts always emitted |

## Build scaling by BAM shape

Selected cases; build at 8 threads. `Mreads/s` is accepted records per second. `Features` is the sum of `units[].total_features`; for `both` and group layouts this is a logical count, not the number of stored spans. `Load` is a `--list-units` process and `Query` is 5,000 BED windows including process startup and load.

| Case | Reads | Contigs | Contig len | Read len | CIGAR | BAM MiB | Features | Build s | Mreads/s | Index MiB | Runs | Build RSS MiB | Load ms | Query ms |
|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| base_800k | 800,000 | 4 | 10,000,000 | 150 | m | 41.9 | 800,000 | 0.374 | 2.14 | 2.71 | 39,154 | 11.9 | 7.4 | 20.1 |
| contigs_1 | 800,000 | 1 | 10,000,000 | 150 | m | 41.4 | 800,000 | 1.247 | 0.64 | 1.88 | 6 | 15.1 | 6.9 | 18.7 |
| contigs_64 | 800,000 | 64 | 1,000,000 | 150 | m | 40.9 | 800,000 | 0.237 | 3.38 | 4.29 | 121,857 | 12.2 | 10.4 | 31.2 |
| contigs_256 | 800,000 | 256 | 1,000,000 | 150 | m | 41.0 | 800,000 | 0.300 | 2.66 | 11.98 | 499,297 | 18.4 | 11.8 | 29.0 |
| contigs_4096 | 4,096 | 4,096 | 100,000 | 150 | m | 0.2 | 4,096 | 0.670 | 0.01 | 1.41 | 4,096 | 21.9 | 43.9 | 85.4 |
| contigs_8192 | 8,192 | 8,192 | 100,000 | 150 | m | 0.5 | 8,192 | 2.049 | 0.00 | 2.82 | 8,192 | 38.3 | 85.5 | 154.1 |
| length_500mb | 800,000 | 4 | 500,000,000 | 150 | m | 44.0 | 800,000 | 0.531 | 1.51 | 16.69 | 753,552 | 39.5 | 6.6 | 22.5 |
| length_max_2gb | 50,000 | 1 | 2,147,483,647 | 150 | m | 2.9 | 50,000 | 0.143 | 0.35 | 1.16 | 49,812 | 12.5 | 6.3 | 18.9 |
| reads_10k | 10,000 | 4 | 10,000,000 | 150 | m | 0.5 | 10,000 | 0.022 | 0.45 | 0.21 | 9,616 | 12.4 | 8.2 | 18.8 |
| reads_4m | 4,000,000 | 4 | 10,000,000 | 150 | m | 206.6 | 4,000,000 | 1.677 | 2.39 | 9.40 | 5 | 18.2 | 9.9 | 20.0 |
| readlen_50 | 800,000 | 4 | 10,000,000 | 50 | m | 19.8 | 800,000 | 0.333 | 2.41 | 7.04 | 291,560 | 12.5 | 6.9 | 29.4 |
| readlen_1000 | 800,000 | 4 | 10,000,000 | 1,000 | m | 217.0 | 800,000 | 0.945 | 0.85 | 2.00 | 4 | 12.8 | 6.4 | 19.7 |
| readlen_10000 | 80,000 | 4 | 10,000,000 | 10,000 | m | 201.7 | 80,000 | 0.820 | 0.10 | 0.31 | 4 | 12.5 | 7.8 | 18.8 |
| cigar_skips | 800,000 | 4 | 10,000,000 | 150 | skips | 42.1 | 800,000 | 0.420 | 1.91 | 3.45 | 78,188 | 12.5 | 5.5 | 26.8 |
| cigar_mixed | 800,000 | 4 | 10,000,000 | 150 | mixed | 48.2 | 800,000 | 0.520 | 1.54 | 2.84 | 45,642 | 12.5 | 6.4 | 20.2 |
| cigar_manyops | 40,000 | 4 | 10,000,000 | 1,000 | manyops | 7.5 | 40,000 | 0.402 | 0.10 | 0.57 | 22,004 | 20.0 | 7.3 | 20.4 |
| samples_8 | 800,000 | 4 | 10,000,000 | 150 | m | 42.9 | 1,600,000 | 0.488 | 1.64 | 13.26 | 588,346 | 20.9 | 8.7 | 21.1 |
| samples_64 | 250,000 | 4 | 10,000,000 | 150 | m | 13.6 | 250,000 | 0.345 | 0.72 | 5.58 | 246,418 | 14.4 | 12.9 | 25.3 |
| samples_256 | 100,000 | 2 | 10,000,000 | 150 | m | 5.5 | 100,000 | 0.758 | 0.13 | 2.60 | 99,714 | 13.0 | 26.5 | 36.2 |
| unknown_rg | 800,000 | 4 | 10,000,000 | 150 | m | 43.0 | 800,000 | 0.382 | 2.09 | 10.76 | 435,210 | 20.0 | 10.7 | 28.6 |
| no_rg | 800,000 | 4 | 10,000,000 | 150 | m | 41.7 | 800,000 | 0.360 | 2.22 | 2.73 | 39,308 | 20.0 | 7.6 | 20.1 |
| clustered_90 | 800,000 | 4 | 10,000,000 | 150 | m | 40.8 | 800,000 | 0.360 | 2.22 | 3.01 | 56,783 | 20.0 | 6.0 | 19.6 |
| unmapped_tail | 900,000 | 4 | 10,000,000 | 150 | m | 46.2 | 800,000 | 0.394 | 2.03 | 2.71 | 39,154 | 21.7 | 6.9 | 21.2 |
| single_contig_many_runs | 800,000 | 1 | 500,000,000 | 150 | m | 42.7 | 800,000 | 1.404 | 0.57 | 14.58 | 628,689 | 51.3 | 7.2 | 22.6 |

### Observations

- **Format 1.1.0 has two payloads.** Coverage grows with canonical runs and is unchanged from 1.0.0. Feature counts grow with accepted alignments. On the 800,000-read cases, feature spans add about 2 MiB to the index and up to about 0.1 s to the build; the coverage run tables and `Load`/`Query` times are essentially unchanged.
- **Coverage saturation no longer collapses the whole index.** Four million 150 bp reads over 40 Mb produce only five coverage runs, but they still store 4,000,000 features and a 9.40 MiB index. In 1.0.0 the same shape was a near-empty index. Feature storage is the price of column 4.
- **Index size is about 19 bytes per run plus 2.5-3.2 bytes per stored feature and 0.26 KiB per contig.** `cigar_skips` (`N` splits blocks) produces roughly twice the runs and a larger coverage payload, while its feature payload stays at one span per read. Sparse long-read data stores fewer features and many fewer runs, so its index is smaller than the read count alone would suggest.
- **Contigs are the task unit.** The same 800,000 reads take 0.37 s across four 10 Mb contigs and 1.25 s on one 10 Mb contig. One contig cannot use more than one worker; `--threads` does not help a single-contig BAM.
- **Read length changes the BAM/CPU trade-off.** 50 bp reads fragment into 292k runs and a 7.04 MiB index; 1,000 and 10,000 bp reads saturate the test contigs at four runs, but their 800k and 80k features keep the index at 2.00 and 0.31 MiB. Build time is dominated by BAM bytes.
- **Many samples and units scale index size with stored features.** Eight samples with `--layout both` produce 9 units, 588k aggregate runs, 800k stored feature spans (1.6M logical), and a 13.26 MiB index at 20.9 MiB RSS. 256 samples produce 256 unit files and a 2.60 MiB index at 13.0 MiB RSS. Build time includes unit-file finalization and grows with unit count.
- **Unknown or absent read groups change unit structure.** Declaring a separate `--unknown-sample` splits reads into another unit, increasing aggregate runs and index size; assigning all orphans to one unit is cheaper. Feature counts are unchanged because every accepted read still contributes one span.
- **Unmapped tails are free.** 100,000 unmapped records increased BAM size but not examined/accepted counts, features, runs, or index size.
- **CSI at `INT32_MAX` works.** A 2,147,483,647 bp contig built and queried correctly with 50,000 features; a 500 Mb contig used BAI, and explicit CSI at 10 Mb was byte-identical in payload to BAI.

## Thread scaling (847 MiB cohort, 20,456,242 alignments, `--layout both`)

| Threads | Build s | Speedup | Efficiency | Build RSS MiB | Reads/s |
|---:|---:|---:|---:|---:|---:|
| 1 | 39.81 | 1.00x | 100% | 97.2 | 0.51M |
| 2 | 20.58 | 1.93x | 97% | 109.3 | 0.99M |
| 4 | 11.82 | 3.37x | 84% | 196.5 | 1.73M |
| 8 | 7.17 | 5.56x | 69% | 371.4 | 2.85M |
| 16 | 5.14 | 7.74x | 48% | 720.7 | 3.98M |
| 30 | 4.34 | 9.17x | 31% | 1,319.9 | 4.71M |

The index is 188.03 MiB at every thread count and all file hashes match. Peak RSS grows by roughly 42 MiB per worker, from 97 MiB at one thread to 1.32 GiB at 30.

Feature spooling and sorting still add cost: compared with the coverage-only 0.1.1 run on the same BAM, one thread rises from 37.77 s to 39.81 s, eight threads from 6.57 s to 7.17 s, and 30 threads from 3.91 s to 4.34 s. That is about 5-11% more build time for 20.5M stored feature spans (60.7 MiB of feature data), down from 12-17% before the optimization pass below. Scaling efficiency is unchanged because feature work is parallelized per contig.

## Optimization pass

After the first format 1.1.0 measurements, the build and query paths were profiled and changed. These numbers compare the unoptimized 1.1.0 binary with the optimized one on the same hardware, BAMs, and 188.03 MiB index. Query timings here write to `/dev/null`; the stress-matrix tables above include output files and are slightly slower.

| Operation | Before | After | Change |
|---|---:|---:|---:|
| Build 847 MiB cohort, 8 threads, tmpfs | 7.38 s | 7.17 s | -3% |
| Build 847 MiB cohort, 8 threads, ext4 | 8.5-8.8 s | 7.2-7.3 s | about -15% |
| Build 800k-read synthetic BAM | 0.414 s | 0.374 s | -10% |
| BED, 10,789 rows, with feature counts | 0.118 s | 0.081 s | -32% |
| BED, 2,789 rows, with feature counts | 0.054 s | 0.043 s | -20% |
| BED, 100,789 rows, with feature counts | 0.78 s | 0.44 s | -44% |
| BED, 500,789 rows, coverage-only path | 1.54 s | 1.15 s | -26% |
| BED, 500,789 rows, with feature counts | 3.79 s | 1.48 s | -61% |
| 1 kb tiling, 603,819 rows, with feature counts | 4.1-4.5 s | 1.20 s | about -72% |

Changes:

- **One count pass per BED row.** BED records are validated before output is opened, and the coverage and feature counts from that pass are now retained for `write_row` instead of being recomputed.
- **Allocation-free index lookups.** Name maps use heterogeneous `string_view` lookup, and `count()`/`feature_count()` resolve the unit and contig once. Derived units precompute member sample IDs, so a pooled query no longer allocates and hashes five member names per row.
- **Block-first feature fast path.** Before decoding a feature block, the reader checks its first value. Windows in a gap before the block return their rank immediately. Little-endian hosts also binary-search the block-last array directly.
- **Batch feature counts.** `CovIndex::feature_counts()` decodes each feature array once and answers many windows on one contig. `covidx-query` uses it for BED and tiling inputs with at least 250,000 rows; smaller inputs keep the cheaper per-window block lookup.
- **Buffered spool reads.** Run and feature spool readers now read 64 KiB at a time instead of calling `fread` for every 8- or 16-byte record.
- **No fsync on temporary spools.** Spools are process-local intermediates that are read immediately after close; unit files and manifests are still fsynced and CRC-checked when published.
- **Skip a redundant sort.** Feature ends are sorted only when they are not already non-decreasing, which is common for fixed-length reads.

Every optimized build produced byte-identical index files when `SOURCE_DATE_EPOCH=0` was set. Query output was byte-identical on the 2,789, 500,789, and 603,819-row cases. The batched feature path raises peak RSS from 216 to 259 MiB for 500k BED rows because it holds two decoded feature arrays and the result vector; the 250,000-row threshold keeps small queries on the lower-memory path.

## Earlier scaling fixes

The fixes below were made before feature-count storage was added. The measurements in this section are historical; the current results appear above. The code paths are unchanged.

### 1. Per-contig BAM header/index reloads

Every contig task constructed a `BamRegionReader`, which re-read the full BAM header and re-loaded the full index. On an 8,192-contig BAM the header took 1.8 ms and the index 10.7 ms per contig; iterating a contig with shared metadata took 0.14 ms. The reader now caches immutable headers/indexes per worker thread while each reader keeps its own file handle and iterator.

| Case | Before | After | Change |
|---|---:|---:|---:|
| 8,192 contigs, 1 thread | 93.34 s | 3.34 s | -96% |
| 8,192 contigs, 8 threads | 13.51 s | 1.67 s | -88% |
| 8,192 contigs, 30 threads | 5.23 s | 1.78 s | -66% |
| 4,096 contigs, 8 threads | 3.52 s | 0.53 s | -85% |
| 20.5M-read cohort, 1 thread | 95.36 s | 37.77 s | -60% |
| 20.5M-read cohort, 30 threads | 4.83 s | 3.91 s | -19% |

Output from the new reader was byte-identical to the old index at 1, 8, and 30 threads.

### 2. Quadratic manifest scope validation

`parse_manifest` linearly searched all contigs for every scope interval. Since a full build has one scope interval per contig, loading was O(contigs²). It now builds a name→ID map once.

| Index | Before | After | Change |
|---|---:|---:|---:|
| 4,096 contigs | 0.090 s | 0.041 s | -54% |
| 8,192 contigs | 0.301 s | 0.083 s | -72% |

The current format 1.1.0 loads are 0.044 s for 4,096 contigs and 0.086 s for 8,192 contigs, still approximately linear.

## Build policy effects (175 MiB, 3,534,590-read BAM, 8 threads)

| Variant | Build s | Index MiB | Runs | Features | Covered bp | Accepted |
|---|---:|---:|---:|---:|---:|---:|
| pooled_default | 1.484 | 26.30 | 875,555 | 3,534,590 | 158,211,560 | 3,534,590 |
| per_sample | 1.633 | 26.30 | 875,555 | 3,534,590 | 158,211,560 | 3,534,590 |
| both | 1.719 | 41.86 | 1,751,110 | 7,069,180 | 316,423,120 | 3,534,590 |
| min_mapq_20 | 1.486 | 26.30 | 875,555 | 3,534,590 | 158,211,560 | 3,534,590 |
| min_mapq_30 | 1.388 | 21.53 | 716,079 | 2,793,272 | 128,712,423 | 2,793,272 |
| min_mapq_60 | 0.966 | 0.06 | 0 | 0 | 0 | 0 |
| min_mapq_255 | 0.972 | 0.06 | 0 | 0 | 0 | 0 |
| flags_none | 1.508 | 26.30 | 875,555 | 3,534,590 | 158,211,560 | 3,534,590 |
| flags_plus_dups | 1.499 | 26.30 | 875,555 | 3,534,590 | 158,211,560 | 3,534,590 |
| include_supplementary | 1.508 | 26.30 | 875,555 | 3,534,590 | 158,211,560 | 3,534,590 |
| no_deletions | 1.521 | 26.64 | 885,920 | 3,534,590 | 158,188,329 | 3,534,590 |
| ambiguous_exclude | 1.716 | 26.30 | 875,555 | 3,534,590 | 158,211,560 | 3,534,590 |
| checksum_source | 3.256 | 26.30 | 875,555 | 3,534,590 | 158,211,560 | 3,534,590 |
| contigs_one | 0.533 | 0.97 | 32,627 | 132,597 | 5,869,292 | 132,597 |
| scope_random | 0.322 | 0.01 | 12 | 89 | 2,428 | 89 |

- This BAM has no duplicates, secondary, supplementary, or QC-fail records, so flag variants are intentionally identical. Its maximum MAPQ is 42, which is why MAPQ 60 and 255 accept nothing while still scanning every record.
- MAPQ 30 removes 21% of reads, 19% of covered bases, 21% of features, and 18% of index size.
- `--layout both` doubles the coverage payload but stores feature spans once in the sample units; the pooled unit derives its feature counts. The total index grows from 26.30 to 41.86 MiB, not to 52.60 MiB.
- `--no-count-deletions` splits at `D`, adding 1.2% runs and 0.34 MiB. Feature spans still use the full `M/D/N/=/X` reference span, so the feature count is unchanged.
- `--ambiguous-read-bases exclude` costs 16% build time on this BAM even though it contains no `N` bases, because every aligned read base is inspected. Feature spans are unchanged.
- `--checksum-source` adds one full BAM pass: 1.77 s here.
- `--regions` is cheaper than `--contigs` when only a small part of a contig is needed; the 80-interval scope examined 89 reads and built in 0.32 s.

## Query behavior

Index used: 188.03 MiB, 7,075,359 runs, 6 units, 20,456,242 stored feature spans. `Query` is 5,000 BED windows plus process startup; "load" is a separate `--list-units` process.

| Load mode | Load s | Load RSS MiB | Query s | Query RSS MiB |
|---|---:|---:|---:|---:|
| mmap + CRC + structural | 0.313 | 80.9 | 0.348 | 150.5 |
| mmap, no CRC, structural | 0.023 | 80.8 | 0.056 | 150.5 |
| heap, no CRC, structural | 0.894 | 251.3 | 1.014 | 250.7 |
| mmap, no CRC, validation none | 0.018 | 25.4 | 0.051 | 110.5 |
| mmap, no CRC, validation full | 0.697 | 193.0 | 0.692 | 194.2 |

- mmap avoids both the 0.9 s heap copy and about 170 MiB of RSS.
- CRC costs about 0.29 s on 188 MiB. Full validation costs about 0.65 s because it reads every coverage run and every feature block.
- Structural validation has no measurable cost over `none` for this workload (0.023 vs 0.018 s load).
- Query cost is output-bound for many small windows: 1.8k rows 0.034 s, 10.8k 0.088 s, 100.8k 0.467 s, 500.8k 1.898 s with feature counts (about 3.8 µs/row at scale). RSS rose to 259 MiB for 500k BED rows because the BED is materialized and the batch feature path holds decoded arrays before queries run.
- Tiling: 604k 1 kb windows on one contig took 1.491 s with feature counts (2.5 µs/row); 21k 1 Mb windows across all contigs took 0.256 s; 10 Mb and 100 Mb tilings are startup-dominated (0.060 s and 0.033 s for 2,369 and 462 rows).
- Run density: 875k runs and 3.5M features in a 26.30 MiB index loaded in 0.012 s and queried 5.8k rows in 0.037 s; 7.08M runs and 20.5M stored features in 188.03 MiB loaded in 0.032 s and queried in 0.064 s. An 8x larger run table plus 20x more feature data costs about 1.7x per query thanks to the bucket tables, feature block search, and lazy mmap.

### Feature queries

Feature counts are part of every query. A 2,789-window mixed-width BED takes about 0.04 s on the per-window path (about 2.5 µs per row). The feature path locates and decodes one block per feature array, so the cost is per query rather than per covered base. The larger working set and page faults into the feature sections contribute on a cold mmap.

The pooled unit was compared with the sum of its five sample units over the same 2,789 windows: 0 mismatches. This exercises the derived-unit path that group and `both`-layout pooled units use, where feature counts are summed from member samples without duplicate storage.

A full-validation load of the six-unit index took 0.712 s with 193 MiB RSS, including every feature block and span/order check.

For large inputs, `covidx-query` batches feature counts. A 500,789-row BED took 1.936 s including output, and a 603,819-row 1 kb tiling took 1.491 s, versus 4.1-4.5 s for per-window lookups before the batch path. The batch path keeps peak RSS at 256 MiB for the 500k BED case.

## Format 1.1.0 index summary

`bench/INDEX_SUMMARY.md` has the full unit-by-unit breakdown. The two reference builds are:

| Index | BAM MiB | Runs | Stored features | Index MiB | Coverage MiB | Feature MiB | B/run | B/feature |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| five-sample, `both` | 847.4 | 7,075,359 | 20,456,242 | 188.03 | 127.17 | 60.73 | 18.85 | 3.11 |
| single-sample, `pooled` | 174.8 | 875,555 | 3,534,590 | 26.30 | 15.55 | 10.69 | 18.62 | 3.17 |

The coverage portion of the five-sample index matches the 127.3 MiB measured for the same BAM under format 1.0.0. Feature storage accounts for the entire increase. The pooled unit in the `both` build reports 20,456,242 logical features but stores no feature bytes; summing `units[].total_features` is therefore not a byte estimate.

## Correctness checks

- 13 independent `samtools` + plain-Python CIGAR reference cases passed: `M`, `D`, `N`, `I`, `S`, `H`, `=`/`X`, mixed, 600-op, deletion-as-gap, ambiguous-base cover/exclude, MAPQ/flag filtering, supplementary inclusion. Each compared 253 windows for both column 5 and column 4, including 10,000 feature spans on `chr1` for the unfiltered cases, 5,300 after MAPQ/flag filtering, and 7,900 with supplementary alignments.
- 52 independent indexes built across all shapes passed `--validation full` and emitted 0 inconsistent coverage rows and 0 inconsistent feature rows. Every index query emitted the feature column and was checked for non-negative counts.
- 10 determinism checks passed: five synthetic shapes at 1 vs 8 threads, and the 20.5M-read cohort at 1 vs 2, 4, 8, 16, and 30 threads. All index files were byte-identical, including feature segments.
- The optimized build produced byte-identical index files with `SOURCE_DATE_EPOCH=0`, and the optimized query produced byte-identical output on 2,789, 500,789, and 603,819-row cases.
- The unit suite has 176 test cases, including batch feature-count coverage; the full suite passes under GCC Debug with `-Werror`, Clang 19, ASan+UBSan, and TSan.
- The five libFuzzer targets cover BED, CIGAR, JSON, manifest, and run-index query parsing. They are not part of this matrix run; feature segment corruption is covered by the unit and corruption tests.

## Remaining limits

- **Feature storage grows with accepted alignments.** A saturated region can have very few coverage runs and still store millions of feature spans. The four-million-read case is 9.40 MiB even though its coverage table has five runs.
- **Feature counts add query cost.** Small queries pay about 2.5 µs per row; large queries use the batch path at the cost of decoding each feature array once and about 40 MiB extra RSS on a 500k-row BED.
- **Full validation is expensive on feature-heavy indexes.** It reads every feature block and span. Use structural validation for routine loads and reserve `--validation full` for diagnostics or after copying an index.
- **Single-contig BAMs do not parallelize.** Contig-level tasks mean a one-contig, 800k-read BAM takes 1.25 s while the same reads across four contigs take 0.37 s.
- **Many-contig finalization is partly serial.** Unit finalization runs per unit; a pooled-only build with thousands of contigs is one unit, so 8,192 contigs still take 2.05 s at 8 threads.
- **Per-contig index overhead is large for tiny contigs.** 4,096/8,192-contig indexes use about 0.26 KiB per contig; a 100 kb contig with one read is much more expensive than the read data itself.
- **BED queries materialize all rows.** 500k rows consume about 216 MiB RSS before the batch feature pass adds its working set; chunk very large BEDs if memory is constrained.
- **Heap loading is expensive.** Prefer mmap and `--no-validate-crc` for trusted, repeatedly opened indexes; reserve `--validation full` for diagnostics.
- **30-worker builds need about 1.3 GiB RSS** for an 847 MiB BAM, roughly 42 MiB per worker plus baseline.
