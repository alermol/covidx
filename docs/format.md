# covidx on-disk format, version 1.1.0

This document is normative. The implementation constants live in `include/covidx/format.h`.

Format `1.1.0` is the only format read or written by covidx 0.2.0. Every index stores feature-count segments in addition to coverage runs. Indexes with any other format version are rejected.

## Directory layout

```
<name>.covidx/
  manifest.json
  contigs.tsv
  samples.tsv
  units.tsv
  contigs.tsv.crc32
  samples.tsv.crc32
  units.tsv.crc32
  data/
    u000000.bin
    u000000.bin.crc32
    u000001.bin
    u000001.bin.crc32
    ...
```

Unit files are named by zero-padded numeric unit ID, never by user-controlled names. One unit file contains all contigs for that unit.

## Constants

| Constant | Value |
|---|---|
| Format version | `1.1.0` |
| Unit magic | `COVIDXU1` (8 bytes, no terminator) |
| Unit header size | 64 bytes |
| Metadata record size | 48 bytes |
| Pack alignment | 8 bytes |
| Contig length limit | `INT32_MAX` (2,147,483,647) |
| Byte order | little-endian |
| Feature segment magic | `COVIDXF1` (8 bytes, no terminator) |
| Feature segment header size | 64 bytes |
| Feature block size | 64 values |
| Unit header flag | bit 0: feature-count segments present |

Positions are always `uint32`; `cum` entries are `uint64`; bucket-table entries are `uint32`. Feature block offsets are `uint64`; feature values and block-last entries are `uint32`; feature values use unsigned LEB128 varints.

## Unit file header (64 bytes)

| Offset | Type | Field |
|---:|---|---|
| 0 | `char[8]` | magic `COVIDXU1` |
| 8 | `u32` | unit ID |
| 12 | `u32` | contig count `N` |
| 16 | `u32` | metadata record size (48) |
| 20 | `u32` | flags: bit 0 marks feature-count segments |
| 24 | `u64` | metadata offset (64, 8-aligned) |
| 32 | `u64` | data offset = `align8(64 + N*48)` |
| 40 | `u64` | reserved (0) |
| 48 | `u64` | reserved (0) |
| 56 | `u64` | reserved (0) |

Metadata records follow, one per indexed contig, in contig-ID order.

## Metadata record (48 bytes)

| Offset | Type | Field |
|---:|---|---|
| 0 | `u64` | segment file offset; 0 if `R=0` |
| 8 | `u64` | run count `R` |
| 16 | `u64` | total covered positions inside the scope |
| 24 | `u64` | maximum bucket population |
| 32 | `u32` | bucket width exponent `w` (6..20; 20 if `R=0`) |
| 36 | `u8` | flags (0) |
| 37 | `u8` | reserved (0) |
| 38 | `u16` | reserved (0) |
| 40 | `u64` | feature segment offset; 0 if this contig has no features |

## Segment layout

Segments start at an 8-byte-aligned offset. Arrays are written in this order:

```
start[R]    uint32, 4 bytes each
end[R]      uint32, 4 bytes each
cum[R+1]    uint64, 8 bytes each  (R+1 entries)
first[B]    uint32, 4 bytes each
pad to 8
```

`B = ceil(L / 2^w) + 1` for `R > 0`; `B = 0` for `R = 0`.

The bucket table satisfies `first[b] = number of runs with start < b*2^w`. Prefix coverage is:

```
C(x) = 0                                         if x == 0, R == 0, or no start <= x
C(x) = cum[j] + min(x, end[j]) - start[j]        otherwise, j = last run with start[j] <= x
C(x) = total_covered                             if x >= L
covered_count(X,Y) = C(Y) - C(X)
```

`w` is the smallest value in `[6,20]` with `floor(L / 2^w) <= R` (clamped to 20). Runs are canonical: sorted, `start < end`, and each run ends strictly before the next run starts. The bucket search therefore examines at most `2^19` starts (19 comparisons).

When `R = 0`, `segment_offset = 0`, `total_covered = 0`, `max_bucket_population = 0`, and `w = 20`; no coverage segment is present. A contig can still have a feature segment when `R = 0`.

## Feature segment

Every format `1.1.0` build sets unit-header flag bit 0 on sample units and on a pooled-only unit. Each contig with at least one accepted alignment has a feature segment at the 8-byte-aligned `feature_offset`. The segment contains two independent sorted arrays: feature starts and feature ends. A feature is one accepted alignment over its full reference span `[POS, POS + M/D/N/=/X)`. Starts and ends are sorted separately, so the overlap count is:

```
feature_count(X,Y) =
    lower_bound(starts, Y) - upper_bound(ends, X)
```

Zero-reference-span alignments have no feature. Group units never store feature segments; their counts are the sum of member sample units. The same applies to the pooled unit in a `both` layout. A pooled-only layout stores features in its single pooled unit.

Segment header (64 bytes):

| Offset | Type | Field |
|---:|---|---|
| 0 | `char[8]` | magic `COVIDXF1` |
| 8 | `u32` | feature segment version (1) |
| 12 | `u32` | block size (64) |
| 16 | `u64` | feature count `F` |
| 24 | `u64` | starts section offset (64) |
| 32 | `u64` | ends section offset |
| 40 | `u64` | starts section bytes |
| 48 | `u64` | ends section bytes |
| 56 | `u64` | reserved (0) |

Each section is:

```
block_offset[B]  u64, 8 bytes each   offset from the feature segment start
block_last[B]    u32, 4 bytes each   last value in each block
block data       varint values
```

`B = ceil(F / 64)`. Block 0 data begins at `section_offset + B*12`. Each block starts with the absolute first value followed by deltas between consecutive values, all encoded as unsigned LEB128. Values in each block are non-decreasing. For the last block, `block_last` is the final value; for earlier blocks it is the final value of that block. The `block_last` array lets a reader locate the first block that can contain a query value with a binary search, then decode only that block. This also handles duplicate values across block boundaries.

When `F = 0` for a contig, `feature_offset = 0`; the unit-header flag still records that the unit stores feature counts.

## Size model

The format has two payloads with different growth rules. Coverage grows with canonical runs, not sequencing depth. Feature counts grow with accepted alignments, because each alignment contributes one start and one end.

For a contig with `R > 0` runs and `B` bucket entries, the coverage segment is:

```
4R + 4R + 8(R + 1) + 4B = 16R + 8 + 4B bytes
```

`B = ceil(L / 2^w) + 1`, and `w` is chosen so `B` is close to `R + 1`. Dense data therefore costs roughly 20 bytes per run.

A feature segment with `F > 0` features contains two sections. Each section has `ceil(F / 64)` block offsets (8 bytes), block-last values (4 bytes), and one absolute value plus sorted deltas as LEB128 varints. The 64-byte segment header and 8-byte alignment padding are additional. Sorted coordinate deltas are usually one or two bytes for short reads, so the observed cost is about 3 bytes per stored feature; sparse or long-read spans can need more.

The unit file also has a 64-byte header, 48 bytes of metadata per contig, and up to 7 bytes of alignment padding per segment.

The five-sample `data-raw/multi-sample.bam` (`--layout both`, 263 contigs, 7,075,359 runs) stores 20,456,242 features in 197,103,360 unit bytes: 133,346,372 coverage bytes (18.85 per run), 63,679,776 feature bytes (3.11 per stored feature), and 77,212 bytes of headers, metadata, and padding. The pooled unit reports `total_features = 20,456,242` in the manifest but stores no feature segment because it is derived from the five sample units. Summing `units[].total_features` therefore counts logical features, not stored bytes. See `bench/INDEX_SUMMARY.md` and `bench/index_summary.py` for measured breakdowns.

## Manifest

`manifest.json` is canonical JSON:

- object keys sorted by unsigned UTF-8 bytes;
- two-space indentation;
- LF line endings;
- exactly one trailing newline;
- integers in decimal without leading zeros;
- no floating-point numbers.

Required top-level keys:

`format`, `format_version`, `builder_version`, `created`, `endianness`, `source`, `policy`, `layout`, `contigs`, `scope`, `samples`, `units`, `files`, `index`, `semantics`.

`format` is `covidx`; `endianness` is `little`; `format_version` must be `1.1.0`. Readers ignore unknown keys but reject any other format version. Duplicate keys are rejected.

`scope` is a sorted array of merged, non-touching intervals, one or more per indexed contig. `index.partial_scope` is true when any indexed contig is not indexed over its full `[0,L)` length. Query windows must be fully contained in one scope interval.

Every `1.1.0` manifest has `units[].total_features` for every unit and a non-empty `semantics.features`. For sample and pooled-only units, `total_features` is the stored feature count; for group and `both`-layout pooled units it is the sum of member sample counts.

`files` maps relative POSIX paths to `{crc32, size}` for every unit file and the three TSVs; it does not include sidecars or the manifest itself.

`created` comes from `SOURCE_DATE_EPOCH` when set, otherwise the current UTC time. No mtime, inode, hostname, thread count, or scheduling data is stored.

## TSVs

`contigs.tsv`: `<id>\t<name>\t<length>\n`, indexed contigs in ID order.

`samples.tsv`: `<sample>\t<rg_id>\n`, sorted by sample then RG ID. A sample with no read groups has one row with `-`.

`units.tsv`: `<unit_id>\t<kind>\t<unit_name>\t<member>\n`, sorted by unit ID then member. `kind` is `sample`, `group`, or `pooled`; pooled uses `-` as the member.

The manifest is authoritative. Loaders compare the TSVs against it and reject mismatches.

## Checksums

CRC-32/ISO-HDLC is the zlib `crc32` algorithm: polynomial `0xEDB88320`, initial and final xor `0xFFFFFFFF`. The test vector `"123456789"` gives `0xCBF43926`.

Each unit file and TSV has a sidecar `<path>.crc32` containing exactly eight lowercase hexadecimal digits and one LF. The manifest CRC is authoritative. A sidecar that disagrees with the manifest is a `FormatError`.

## Save and load

`save()` writes a temporary directory, copies or hardlinks immutable unit files, fsyncs files and directories, verifies sizes and CRCs, and publishes with `rename`. Existing indexes are moved aside and restored if publishing fails.

`load()` validates the manifest, file sizes, CRCs (unless disabled), sidecars, TSVs, unit headers, metadata, coverage arrays, and feature segments according to the requested validation level:

- `none`: headers, offsets, sizes, and manifest/TSV consistency;
- `structural` (default): sampled array invariants, including the first and last feature blocks;
- `full`: every run, `cum`, bucket entry, feature block, and feature span/order check.

## Big-endian hosts

The on-disk format is always little-endian. On big-endian hosts the reader copies and byte-swaps arrays; correctness is required, performance is secondary.
