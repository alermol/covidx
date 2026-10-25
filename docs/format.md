# covidx on-disk format, version 1.0.0

This document is normative. The implementation constants live in `include/covidx/format.h`.

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
| Format version | `1.0.0` |
| Unit magic | `COVIDXU1` (8 bytes, no terminator) |
| Unit header size | 64 bytes |
| Metadata record size | 48 bytes |
| Pack alignment | 8 bytes |
| Contig length limit | `INT32_MAX` (2,147,483,647) |
| Byte order | little-endian |

Positions are always `uint32`; `cum` entries are `uint64`; bucket-table entries are `uint32`.

## Unit file header (64 bytes)

| Offset | Type | Field |
|---:|---|---|
| 0 | `char[8]` | magic `COVIDXU1` |
| 8 | `u32` | unit ID |
| 12 | `u32` | contig count `N` |
| 16 | `u32` | metadata record size (48) |
| 20 | `u32` | flags (0) |
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
| 40 | `u32` | reserved (0) |
| 44 | `u32` | reserved (0) |

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

When `R = 0`, `segment_offset = 0`, `total_covered = 0`, `max_bucket_population = 0`, and `w = 20`; no segment is present.

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

`format` is `covidx`; `endianness` is `little`; `format_version` is `MAJOR.MINOR.PATCH`. Readers accept the same major version and ignore unknown keys. Duplicate keys are rejected.

`scope` is a sorted array of merged, non-touching intervals, one or more per indexed contig. `index.partial_scope` is true when any indexed contig is not indexed over its full `[0,L)` length. Query windows must be fully contained in one scope interval.

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

`load()` validates the manifest, file sizes, CRCs (unless disabled), sidecars, TSVs, unit headers, metadata, and arrays according to the requested validation level:

- `none`: headers, offsets, sizes, and manifest/TSV consistency;
- `structural` (default): sampled array invariants;
- `full`: every run, `cum`, and bucket entry.

## Big-endian hosts

The on-disk format is always little-endian. On big-endian hosts the reader copies and byte-swaps arrays; correctness is required, performance is secondary.
