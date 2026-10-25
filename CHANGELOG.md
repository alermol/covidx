# Changelog

## 0.2.0 — 2026-10-25

- Store feature counts in every index, matching default `bedtools coverage` column 4 (`num_features`).
- Store one feature per accepted alignment over its full reference span (`M/D/N/=/X`), independent of deletion and ambiguous-base coverage policy.
- Make `covidx-query` always emit `num_features` before `covered`; there is no opt-in flag.
- Add `CovIndex::feature_count()` and the batch `CovIndex::feature_counts()`.
- Introduce format `1.1.0` with block-varint feature segments as the only supported format. Format `1.0.0` indexes are not readable; rebuild them.
- Group and pooled feature counts sum their member sample counts without duplicating feature storage.
- Extend parity harnesses and tests with column-4 validation.
- Speed up large queries with cached BED counts, allocation-free index lookups, a feature block-first fast path, and batch feature counts for large BED/tiling inputs.
- Speed up build finalization with buffered spool reads, no fsync on temporary spools, and skipping redundant feature-end sorting.
- Document the format `1.1.0` size model, add `bench/index_summary.py`, and refresh the stress matrix for feature counts.
- Remove unused internal accessors, the unreleased `CovIndex::num_features()` alias, a dead JSON parser store, and unused includes; no behavior change.

## 0.1.1 — 2026-10-25

- Make real-BAM parity validation non-vacuous with bounded regions, sample filtering, index policy checks, and explicit empty-run handling.
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
