# Changelog

## 0.1.0 — 2026-10-25

- Initial C++23 coverage index implementation.
- Exact covered-position queries with a canonical run table and bucket table.
- BAM header/index validation, sample demultiplexing, per-sample, pooled, and fixed-group units.
- `--contigs` and `--regions` build scope with strict query containment.
- Deterministic, self-contained unit files, canonical manifest, CRC32 sidecars, atomic save, and validation levels.
- `covidx-build` and `covidx-query` command-line tools.
- Samtools and bedtools parity harnesses, property tests, and golden index.
