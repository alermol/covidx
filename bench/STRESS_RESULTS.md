# covidx stress-test results

Measured 2026-10-24 on Debian trixie, 32 vCPU, 117 GiB RAM, GCC 14.2.0, htslib 1.21, Release build. BAMs, indexes, and temporary work directories were on tmpfs (`/tmp`). All builds used `SOURCE_DATE_EPOCH=0` unless noted.

Reproduce with:

```bash
cmake --build build-release -j
./bench/stress_matrix.py --stage all \
  --work-dir /tmp/covidx-stress \
  --windows 5000 --threads 8
```

The harness writes `stress-results.json` plus one TSV per stage. The results from this run are copied to `bench/stress-results/`.

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
| Query | mmap/heap, CRC on/off, validation none/structural/full, BED 1k-500k rows, 1 kb-100 Mb tiling |

## Build scaling by BAM shape

Selected cases; build at 8 threads. `Mreads/s` is mapped records per second. `Load` is a `--list-units` process (index load) and `Query` is 5,000 BED windows including process startup and load.

| Case | Reads | Contigs | Contig len | Read len | CIGAR | BAM MiB | Build s | Mreads/s | Index MiB | Runs | Build RSS MiB | Load ms | Query ms |
|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| base_800k | 800,000 | 4 | 10,000,000 | 150 | m | 41.9 | 0.313 | 2.56 | 0.73 | 39,288 | 12.2 | 7.5 | 21.7 |
| contigs_1 | 800,000 | 1 | 10,000,000 | 150 | m | 41.4 | 1.149 | 0.70 | 0.00 | 6 | 12.9 | 6.4 | 20.7 |
| contigs_64 | 800,000 | 64 | 1,000,000 | 150 | m | 40.9 | 0.206 | 3.89 | 2.14 | 122,153 | 12.9 | 9.3 | 24.2 |
| contigs_256 | 800,000 | 256 | 1,000,000 | 150 | m | 41.0 | 0.269 | 2.97 | 9.17 | 500,358 | 16.9 | 9.9 | 29.7 |
| contigs_4096 | 4,096 | 4,096 | 100,000 | 150 | m | 0.2 | 0.525 | 0.01 | 1.03 | 4,096 | 21.9 | 41.1 | 90.2 |
| contigs_8192 | 8,192 | 8,192 | 100,000 | 150 | m | 0.5 | 1.673 | 0.00 | 2.07 | 8,192 | 38.0 | 82.9 | 162.8 |
| length_500mb | 800,000 | 4 | 500,000,000 | 150 | m | 44.0 | 0.531 | 1.51 | 13.36 | 753,266 | 39.3 | 7.1 | 26.6 |
| length_max_2gb | 50,000 | 1 | 2,147,483,647 | 150 | m | 2.9 | 0.142 | 0.35 | 0.89 | 49,812 | 13.2 | 7.1 | 22.2 |
| reads_10k | 10,000 | 4 | 10,000,000 | 150 | m | 0.5 | 0.021 | 0.47 | 0.17 | 9,619 | 12.6 | 6.9 | 21.6 |
| reads_4m | 4,000,000 | 4 | 10,000,000 | 150 | m | 206.6 | 1.450 | 2.76 | 0.00 | 6 | 12.7 | 7.7 | 21.3 |
| readlen_50 | 800,000 | 4 | 10,000,000 | 50 | m | 19.8 | 0.282 | 2.84 | 5.04 | 291,277 | 12.7 | 6.3 | 23.5 |
| readlen_1000 | 800,000 | 4 | 10,000,000 | 1,000 | m | 217.0 | 0.894 | 0.90 | 0.00 | 4 | 13.0 | 8.1 | 21.3 |
| readlen_10000 | 80,000 | 4 | 10,000,000 | 10,000 | m | 201.7 | 0.814 | 0.10 | 0.00 | 4 | 12.7 | 6.7 | 20.4 |
| cigar_skips | 800,000 | 4 | 10,000,000 | 150 | skips | 42.1 | 0.388 | 2.06 | 1.46 | 78,198 | 12.7 | 7.3 | 22.7 |
| cigar_mixed | 800,000 | 4 | 10,000,000 | 150 | mixed | 48.2 | 0.464 | 1.73 | 0.84 | 45,341 | 12.7 | 6.7 | 22.6 |
| cigar_manyops | 40,000 | 4 | 10,000,000 | 1,000 | 600-op | 7.5 | 0.397 | 0.10 | 0.41 | 22,004 | 20.9 | 7.2 | 22.4 |
| samples_8 | 800,000 | 4 | 10,000,000 | 150 | m | 42.9 | 0.422 | 1.90 | 10.31 | 588,305 | 22.1 | 7.9 | 23.1 |
| samples_64 | 250,000 | 4 | 10,000,000 | 150 | m | 13.6 | 0.272 | 0.92 | 4.41 | 246,404 | 13.2 | 10.6 | 25.4 |
| samples_256 | 100,000 | 2 | 10,000,000 | 150 | m | 5.5 | 0.714 | 0.14 | 1.99 | 99,714 | 13.2 | 21.1 | 35.2 |
| unknown_rg | 800,000 | 4 | 10,000,000 | 150 | m | 43.0 | 0.350 | 2.29 | 7.99 | 435,206 | 20.9 | 8.7 | 24.5 |
| no_rg | 800,000 | 4 | 10,000,000 | 150 | m | 41.7 | 0.300 | 2.67 | 0.75 | 39,294 | 20.9 | 6.4 | 21.9 |
| clustered_90 | 800,000 | 4 | 10,000,000 | 150 | m | 40.8 | 0.320 | 2.50 | 1.02 | 56,668 | 20.9 | 6.8 | 22.6 |
| unmapped_tail | 900,000 | 4 | 10,000,000 | 150 | m | 46.2 | 0.315 | 2.85 | 0.73 | 39,288 | 21.7 | 7.0 | 22.8 |
| single_contig_many_runs | 800,000 | 1 | 500,000,000 | 150 | m | 42.7 | 1.330 | 0.60 | 11.46 | 628,689 | 37.9 | 7.1 | 25.0 |

### Observations

- **Contigs are the task unit.** The same 800,000 reads take 0.31 s across four 10 Mb contigs and 1.15 s on one 10 Mb contig. One contig cannot use more than one worker; `--threads` does not help a single-contig BAM.
- **Coverage saturation collapses run counts.** Four million 150 bp reads over 40 Mb produce only six runs and a near-empty index: build time tracks reads and BAM size, not runs.
- **Index size is about 19 bytes per run plus about 0.26 KiB per contig.** `cigar_skips` (`N` splits blocks) produces roughly twice the runs and twice the index of plain `M` at the same read count. `no-count-deletions` adds about 1.2% more runs on the real BAM.
- **Read length changes the BAM/CPU trade-off.** 50 bp reads fragment into 291k runs and a 5 MiB index; 1,000 and 10,000 bp reads saturate the test contigs, leaving four runs, while build time is dominated by BAM bytes.
- **Many samples/units scale index size, not peak memory.** 8 samples with `--layout both` produce 9 units, 588k aggregate runs, and a 10.3 MiB index at 22 MiB RSS. 256 samples produce 256 unit files and 2.0 MiB at 13 MiB RSS. Build time includes unit-file finalization and grows with unit count.
- **Unknown/absent read groups change unit structure.** Declaring a separate `--unknown-sample` splits reads into another unit, increasing aggregate runs and index size; assigning all orphans to one unit is cheaper.
- **Unmapped tails are free.** 100,000 unmapped records increased BAM size but not examined/accepted counts or index size.
- **CSI at `INT32_MAX` works.** A 2,147,483,647 bp contig built and queried correctly; a 500 Mb contig used BAI, and explicit CSI at 10 Mb was byte-identical in payload to BAI.

## Thread scaling (847 MiB cohort, 20,456,242 reads, `--layout both`)

| Threads | Build s | Speedup | Efficiency | Build RSS MiB | Reads/s |
|---:|---:|---:|---:|---:|---:|
| 1 | 37.77 | 1.00x | 100% | 78.3 | 0.54M |
| 2 | 19.53 | 1.93x | 97% | 108.2 | 1.05M |
| 4 | 11.15 | 3.39x | 85% | 194.7 | 1.84M |
| 8 | 6.57 | 5.75x | 72% | 367.9 | 3.11M |
| 16 | 4.72 | 8.01x | 50% | 714.0 | 4.34M |
| 30 | 3.91 | 9.65x | 32% | 1,312.6 | 5.23M |

Peak RSS grows by roughly 45 MiB per worker, from 78 MiB at one thread to 1.31 GiB at 30. Scaling is strong through four threads and then flattens as BAM I/O, unit finalization, and serial merge/report phases dominate. The index is byte-identical at every thread count.

## Scaling defects found and fixed

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

Load is now approximately linear (4,096 → 8,192 contigs: 0.041 → 0.083 s).

## Build policy effects (175 MiB, 3,534,590-read BAM, 8 threads)

| Variant | Build s | Index MiB | Runs | Covered bp | Accepted |
|---|---:|---:|---:|---:|---:|
| pooled_default | 1.195 | 15.61 | 875,555 | 158,211,560 | 3,534,590 |
| per_sample | 1.323 | 15.61 | 875,555 | 158,211,560 | 3,534,590 |
| both | 1.451 | 31.18 | 1,751,110 | 316,423,120 | 3,534,590 |
| min_mapq_20 | 1.200 | 15.61 | 875,555 | 158,211,560 | 3,534,590 |
| min_mapq_30 | 1.149 | 12.99 | 716,079 | 128,712,423 | 2,793,272 |
| min_mapq_60 | 0.954 | 0.06 | 0 | 0 | 0 |
| min_mapq_255 | 0.961 | 0.06 | 0 | 0 | 0 |
| flags_none | 1.183 | 15.61 | 875,555 | 158,211,560 | 3,534,590 |
| flags_plus_dups | 1.195 | 15.61 | 875,555 | 158,211,560 | 3,534,590 |
| include_supplementary | 1.218 | 15.61 | 875,555 | 158,211,560 | 3,534,590 |
| no_deletions | 1.197 | 15.95 | 885,920 | 158,188,329 | 3,534,590 |
| ambiguous_exclude | 1.489 | 15.61 | 875,555 | 158,211,560 | 3,534,590 |
| checksum_source | 3.029 | 15.61 | 875,555 | 158,211,560 | 3,534,590 |
| contigs_one | 0.515 | 0.57 | 32,627 | 5,869,292 | 132,597 |
| scope_random | 0.329 | 0.01 | 12 | 2,428 | 89 |

- This BAM has no duplicates, secondary, supplementary, or QC-fail records, so flag variants are intentionally identical. Its maximum MAPQ is 42, which is why MAPQ 60 and 255 accept nothing while still scanning every record.
- MAPQ 30 removes 21% of reads, 19% of covered bases, and 17% of index size.
- `--layout both` exactly doubles the unit payload.
- `--no-count-deletions` splits at `D`, adding 1.2% runs and 2% index size.
- `--ambiguous-read-bases exclude` costs 25% build time on this BAM even though it contains no `N` bases, because every aligned read base is inspected.
- `--checksum-source` adds one full BAM pass (1.8 s here).
- `--regions` is cheaper than `--contigs` when only a small part of a contig is needed; the 80-interval scope examined 89 reads and built in 0.33 s.

## Query behavior

Index used: 127.3 MiB, 7,075,359 runs, 6 units. `Query` is 5,000 BED windows plus process startup; "load" is a separate `--list-units` process.

| Load mode | Load s | Load RSS MiB | Query s | Query RSS MiB |
|---|---:|---:|---:|---:|
| mmap + CRC + structural | 0.210 | 59.9 | 0.234 | 90.7 |
| mmap, no CRC, structural | 0.021 | 59.6 | 0.040 | 90.4 |
| heap, no CRC, structural | 0.616 | 190.2 | 0.612 | 189.8 |
| mmap, no CRC, validation none | 0.009 | 12.5 | 0.033 | 45.4 |
| mmap, no CRC, validation full | 0.251 | 132.4 | 0.269 | 133.2 |

- mmap avoids both the 0.6 s heap copy and about 130 MiB of RSS.
- CRC costs about 0.19 s on 127 MiB; full validation costs about 0.23 s and touches most mapped pages.
- Structural validation has no measurable cost over `none` for this workload (0.040 vs 0.033 s including load).
- Query cost is output-bound for many small windows: 1k rows 0.024 s, 10k 0.056 s, 100k 0.371 s, 500k 1.664 s (about 3.3 µs/row at scale). RSS rose to 187 MiB for 500k BED rows because the BED is materialized before queries run.
- Tiling: 604k 1 kb windows on one contig took 1.08 s (1.8 µs/row); 21k 1 Mb windows across all contigs took 0.065 s; 10 Mb and 100 Mb tilings are startup-dominated (0.034 s and 0.020 s for 2,369 and 462 rows).
- Run density: 875k runs / 15.6 MiB loaded in 0.010 s and queried 5.8k rows in 0.030 s; 7.08M runs / 127.3 MiB loaded in 0.016 s and queried in 0.046 s. An 8x larger run table costs about 1.5x per query thanks to the bucket table and lazy mmap.

## Correctness checks

- 13 independent `samtools` + plain-Python CIGAR reference cases passed: `M`, `D`, `N`, `I`, `S`, `H`, `=`/`X`, mixed, 600-op, deletion-as-gap, ambiguous-base cover/exclude, MAPQ/flag filtering, supplementary inclusion. Each compared 253 windows.
- 52 indexes built across all shapes passed `--validation full` and emitted 0 internally inconsistent rows.
- 15 determinism checks passed: five synthetic shapes at 1 vs 8 threads, and the 20.5M-read cohort at 1 vs 2, 4, 8, 16, and 30 threads. All index files were byte-identical.
- Full unit suite: 154 test cases pass in GCC Debug, GCC `-Werror`, Clang 19, ASan+UBSan, and TSan. All five libFuzzer targets completed clean (94k-452k executions in the final run).

## Remaining limits

- **Single-contig BAMs do not parallelize.** Contig-level tasks mean a one-contig, 800k-read BAM takes 1.15 s while the same reads across four contigs take 0.31 s.
- **Many-contig finalization is partly serial.** Unit finalization runs per unit; a pooled-only build with thousands of contigs is one unit, so 8,192 contigs still take 1.67 s at 8 threads and 3.34 s at one thread.
- **Per-contig index overhead is large for tiny contigs.** 4,096/8,192-contig indexes use about 0.26 KiB per contig; a 100 kb contig with one read is much more expensive than the read data itself.
- **BED queries materialize all rows.** 500k rows consume about 187 MiB RSS; chunk very large BEDs if memory is constrained.
- **Heap loading is expensive.** Prefer mmap and `--no-validate-crc` for trusted, repeatedly opened indexes; reserve `--validation full` for diagnostics.
- **30-worker builds need about 1.3 GiB RSS** for an 847 MiB BAM, roughly 45 MiB per worker plus baseline.
