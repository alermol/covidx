# Changelog

## 0.1.1 — 2026-10-25

- Make real-BAM parity validation non-vacuous with bounded regions, sample
  filtering, index policy checks, and explicit empty-run handling.
- Document read-rich validation data and the limits of the supplied BAMs.
- No library or CLI behavior changes; format remains 1.0.0.

## 0.1.0 — 2026-10-25

- Initial C++23 coverage index implementation.
- Exact covered-position queries with a canonical run table and bucket table.
- BAM header/index validation, sample demultiplexing, per-sample, pooled, and fixed-group units.
- `--contigs` and `--regions` build scope with strict query containment.
- Deterministic, self-contained unit files, canonical manifest, CRC32 sidecars, atomic save, and validation levels.
- `covidx-build` and `covidx-query` command-line tools.
- Samtools and bedtools parity harnesses, property tests, and golden index.
