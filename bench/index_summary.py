#!/usr/bin/env python3
"""Summarize the payload of a covidx 1.1.0 index directory.

Reads manifest.json and the unit headers to split every unit file into
coverage bytes, feature bytes, and fixed overhead. This makes the cost of
feature-count storage visible without relying on the BAM.

Examples:
  ./bench/index_summary.py /tmp/covidx/cohort.covidx
  ./bench/index_summary.py /tmp/covidx/cohort.covidx --json summary.json
"""

import argparse
import json
import struct
from pathlib import Path

UNIT_HEADER = struct.Struct("<8sIIIIQQQQQ")
META_RECORD = struct.Struct("<QQQQIBBHQ")
FEATURE_HEADER = struct.Struct("<8sIIQQQQQQ")


def align8(value):
    return (value + 7) // 8 * 8


def bucket_entries(length, width):
    return (length >> width) + 1


def read_at(handle, offset, size):
    handle.seek(offset)
    data = handle.read(size)
    if len(data) != size:
        raise SystemExit(f"short read at offset {offset}")
    return data


def read_unit(path, contig_lengths):
    with path.open("rb") as handle:
        header = UNIT_HEADER.unpack(
            read_at(handle, 0, UNIT_HEADER.size)
        )
        (
            magic,
            unit_id,
            contig_count,
            meta_size,
            flags,
            metadata_offset,
            data_offset,
            _reserved0,
            _reserved1,
            _reserved2,
        ) = header
        if magic != b"COVIDXU1":
            raise SystemExit(f"{path}: bad unit magic {magic!r}")
        if meta_size != META_RECORD.size:
            raise SystemExit(f"{path}: unsupported metadata record size")
        coverage_bytes = 0
        feature_bytes = 0
        runs_total = 0
        features_total = 0
        contigs_with_features = 0
        for contig_id in range(contig_count):
            offset = metadata_offset + contig_id * meta_size
            record = read_at(handle, offset, meta_size)
            (
                segment_offset,
                runs,
                _total_covered,
                _max_bucket,
                width,
                _flags,
                _reserved8,
                _reserved16,
                feature_offset,
            ) = META_RECORD.unpack(record)
            runs_total += runs
            if runs > 0:
                length = contig_lengths[contig_id]
                coverage_bytes += (
                    runs * 4
                    + runs * 4
                    + (runs + 1) * 8
                    + bucket_entries(length, width) * 4
                )
            if feature_offset == 0:
                continue
            feature = FEATURE_HEADER.unpack(
                read_at(handle, feature_offset, FEATURE_HEADER.size)
            )
            (
                feature_magic,
                version,
                block_size,
                count,
                starts_offset,
                ends_offset,
                starts_bytes,
                ends_bytes,
                _reserved,
            ) = feature
            if feature_magic != b"COVIDXF1":
                raise SystemExit(f"{path}: bad feature magic")
            if version != 1 or block_size != 64:
                raise SystemExit(f"{path}: unsupported feature segment")
            if starts_offset != FEATURE_HEADER.size:
                raise SystemExit(f"{path}: unexpected feature starts offset")
            payload = ends_offset + ends_bytes
            feature_bytes += align8(payload)
            features_total += count
            contigs_with_features += 1
        file_bytes = path.stat().st_size
        overhead = file_bytes - coverage_bytes - feature_bytes
        if overhead < 0:
            raise SystemExit(f"{path}: segment sizes exceed the file size")
    return {
        "file": str(path),
        "bytes": file_bytes,
        "header_bytes": UNIT_HEADER.size,
        "metadata_bytes": align8(
            metadata_offset + contig_count * META_RECORD.size
        ),
        "data_offset": data_offset,
        "flags": flags,
        "runs": runs_total,
        "features": features_total,
        "contigs_with_features": contigs_with_features,
        "coverage_bytes": coverage_bytes,
        "feature_bytes": feature_bytes,
        "overhead_bytes": overhead,
        "feature_bytes_per_feature": round(
            feature_bytes / features_total, 3
        )
        if features_total
        else 0.0,
        "coverage_bytes_per_run": round(coverage_bytes / runs_total, 3)
        if runs_total
        else 0.0,
    }


def summarize(index_dir):
    root = Path(index_dir)
    manifest = json.loads((root / "manifest.json").read_text())
    if manifest.get("format_version") != "1.1.0":
        raise SystemExit(
            f"{root}: expected format 1.1.0, got "
            f"{manifest.get('format_version')!r}"
        )
    contig_lengths = [contig["length"] for contig in manifest["contigs"]]
    units = []
    for unit in manifest["units"]:
        info = read_unit(root / unit["file"], contig_lengths)
        info.update(
            {
                "name": unit["name"],
                "kind": unit["kind"],
                "total_runs": unit["total_runs"],
                "total_features": unit["total_features"],
            }
        )
        units.append(info)
    payload = sum(info["bytes"] for info in units)
    directory = sum(
        path.stat().st_size for path in root.rglob("*") if path.is_file()
    )
    return {
        "index": str(root),
        "format_version": manifest["format_version"],
        "builder_version": manifest["builder_version"],
        "layout": manifest["layout"],
        "contigs": manifest["index"]["contig_count"],
        "units": manifest["index"]["unit_count"],
        "total_runs": manifest["index"]["total_runs"],
        "total_features": sum(
            unit["total_features"] for unit in manifest["units"]
        ),
        "stored_features": sum(info["features"] for info in units),
        "source": manifest["source"],
        "unit_payload_bytes": payload,
        "directory_bytes": directory,
        "coverage_bytes": sum(info["coverage_bytes"] for info in units),
        "feature_bytes": sum(info["feature_bytes"] for info in units),
        "overhead_bytes": sum(info["overhead_bytes"] for info in units),
        "unit_details": units,
    }


def print_summary(summary):
    source = summary["source"]
    print(f"index: {summary['index']}")
    print(
        f"format {summary['format_version']}  "
        f"builder {summary['builder_version']}  "
        f"layout {summary['layout']}"
    )
    print(
        f"source: {source['path']}  "
        f"{source['size'] / (1024 * 1024):.1f} MiB  "
        f"{source['records_accepted']:,} accepted records"
    )
    print(
        f"index: {summary['contigs']:,} contigs  "
        f"{summary['units']:,} units  "
        f"{summary['total_runs']:,} runs  "
        f"{summary['total_features']:,} logical features "
        f"({summary['stored_features']:,} stored)"
    )
    print(
        f"payload: {summary['unit_payload_bytes'] / (1024 * 1024):.2f} MiB "
        f"({summary['directory_bytes'] / (1024 * 1024):.2f} MiB on disk)"
    )
    print(
        f"breakdown: coverage "
        f"{summary['coverage_bytes'] / (1024 * 1024):.2f} MiB, "
        f"features {summary['feature_bytes'] / (1024 * 1024):.2f} MiB, "
        f"overhead {summary['overhead_bytes'] / (1024 * 1024):.2f} MiB"
    )
    print()
    print(
        "unit\tkind\tfile MiB\truns\tfeatures\tcoverage MiB\t"
        "feature MiB\tB/run\tB/feature"
    )
    for unit in summary["unit_details"]:
        print(
            f"{unit['name']}\t{unit['kind']}\t"
            f"{unit['bytes'] / (1024 * 1024):.2f}\t"
            f"{unit['total_runs']:,}\t{unit['total_features']:,}\t"
            f"{unit['coverage_bytes'] / (1024 * 1024):.2f}\t"
            f"{unit['feature_bytes'] / (1024 * 1024):.2f}\t"
            f"{unit['coverage_bytes_per_run']:.2f}\t"
            f"{unit['feature_bytes_per_feature']:.2f}"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("index", help="path to a .covidx directory")
    parser.add_argument("--json", help="write the full summary as JSON")
    args = parser.parse_args()
    summary = summarize(args.index)
    print_summary(summary)
    if args.json:
        Path(args.json).write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
