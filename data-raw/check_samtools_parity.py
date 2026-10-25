#!/usr/bin/env python3
"""Compare covidx covered counts with samtools depth.

The mapping is documented in docs/parity.md. This script deliberately uses
only samtools and the covidx-query CLI, so it also works where pysam is absent.
"""
import argparse
import random
import subprocess
import tempfile
from pathlib import Path


def run(cmd):
    return subprocess.run(cmd, check=True, text=True, capture_output=True)


def contig_length(bam, contig):
    header = run(["samtools", "view", "-H", bam]).stdout
    for line in header.splitlines():
        if not line.startswith("@SQ"):
            continue
        fields = dict(part.split(":", 1) for part in line.split("\t")[1:])
        if fields.get("SN") == contig:
            return int(fields["LN"])
    raise SystemExit(f"contig {contig!r} is not present in {bam}")


def windows_for(length, count):
    windows = []
    for x in range(length):
        for size in range(1, min(50, length - x) + 1):
            windows.append((x, x + size))
    rng = random.Random(17)
    for _ in range(count):
        x = rng.randrange(length)
        windows.append((x, x + 1 + rng.randrange(length - x)))
    return windows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--contig", required=True)
    parser.add_argument("--index", required=True)
    parser.add_argument("--unit", default=None)
    parser.add_argument("--covidx-query", default="covidx-query")
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument(
        "--exclude-flags", type=lambda s: int(s, 0), default=0x704
    )
    parser.add_argument("--include-supplementary", action="store_true")
    parser.add_argument("--no-count-deletions", action="store_true")
    parser.add_argument("--windows-bed", default=None)
    parser.add_argument("--random-windows", type=int, default=500)
    args = parser.parse_args()

    if args.exclude_flags > 0xFFFF:
        raise SystemExit("--exclude-flags does not fit in 16 bits")
    length = contig_length(args.bam, args.contig)

    if args.windows_bed:
        windows = []
        for line in Path(args.windows_bed).read_text().splitlines():
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if fields[0] != args.contig:
                continue
            windows.append((int(fields[1]), int(fields[2])))
    else:
        windows = windows_for(length, args.random_windows)

    with tempfile.TemporaryDirectory(prefix="covidx-samtools-") as tmp:
        bed = Path(tmp) / "windows.bed"
        bed.write_text(
            "".join(f"{args.contig}\t{x}\t{y}\n" for x, y in windows)
        )
        query = [
            args.covidx_query,
            "--index",
            args.index,
            "--bed",
            str(bed),
            "--no-header",
        ]
        if args.unit:
            query += ["--unit", args.unit]
        query_output = run(query).stdout

        effective = (args.exclude_flags | 0x800) & ~(
            0x800 if args.include_supplementary else 0
        )
        default = 0x704
        remove = default & ~effective
        add = effective & ~default
        command = [
            "samtools",
            "depth",
            "-aa",
            "-q",
            "0",
            "-Q",
            str(args.min_mapq),
            "-g",
            str(remove),
            "-G",
            str(add),
        ]
        if not args.no_count_deletions:
            command.append("-J")
        command += ["-r", args.contig, args.bam]
        depth_output = run(command).stdout

    positive = [False] * length
    for line in depth_output.splitlines():
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        position = int(fields[1]) - 1
        positive[position] = int(fields[2]) > 0
    prefix = [0] * (length + 1)
    for i, value in enumerate(positive):
        prefix[i + 1] = prefix[i] + (1 if value else 0)

    rows = [
        line
        for line in query_output.splitlines()
        if line and not line.startswith("#")
    ]
    if len(rows) != len(windows):
        raise SystemExit(
            f"covidx-query returned {len(rows)} rows for {len(windows)} windows"
        )
    for line, (x, y) in zip(rows, windows):
        fields = line.split("\t")
        got = int(fields[4])
        want = prefix[y] - prefix[x]
        if got != want:
            raise SystemExit(
                f"mismatch at {args.contig}:{x}-{y}: "
                f"covidx={got} samtools={want}"
            )
    print(f"ok: {len(windows)} windows on {args.contig}")


if __name__ == "__main__":
    main()
