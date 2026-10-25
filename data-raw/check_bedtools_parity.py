#!/usr/bin/env python3
"""Compare covidx covered counts with default bedtools coverage column 5.

Blocks are generated from samtools depth positive positions, so the covered set
is independent of covidx. See docs/parity.md for the normative procedure.
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

    length = contig_length(args.bam, args.contig)
    rng = random.Random(23)
    if args.windows_bed:
        windows = []
        for line in Path(args.windows_bed).read_text().splitlines():
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if fields[0] == args.contig:
                windows.append((int(fields[1]), int(fields[2])))
    else:
        windows = []
        for x in range(length):
            for size in range(1, min(50, length - x) + 1):
                windows.append((x, x + size))
        for _ in range(args.random_windows):
            x = rng.randrange(length)
            windows.append((x, x + 1 + rng.randrange(length - x)))

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
        if len(fields) >= 3:
            positive[int(fields[1]) - 1] = int(fields[2]) > 0

    with tempfile.TemporaryDirectory(prefix="covidx-bedtools-") as tmp:
        blocks = Path(tmp) / "blocks.bed"
        with blocks.open("w") as out:
            i = 0
            while i < length:
                if not positive[i]:
                    i += 1
                    continue
                start = i
                while i < length and positive[i]:
                    i += 1
                out.write(f"{args.contig}\t{start}\t{i}\n")
        windows_bed = Path(tmp) / "windows.bed"
        windows_bed.write_text(
            "".join(f"{args.contig}\t{x}\t{y}\n" for x, y in windows)
        )
        coverage_output = run(
            ["bedtools", "coverage", "-a", str(windows_bed), "-b", str(blocks)]
        ).stdout
        query = [
            args.covidx_query,
            "--index",
            args.index,
            "--bed",
            str(windows_bed),
            "--no-header",
        ]
        if args.unit:
            query += ["--unit", args.unit]
        query_output = run(query).stdout

    covidx_rows = [
        line
        for line in query_output.splitlines()
        if line and not line.startswith("#")
    ]
    bedtools_rows = [line for line in coverage_output.splitlines() if line]
    if len(covidx_rows) != len(bedtools_rows) or len(covidx_rows) != len(
        windows
    ):
        raise SystemExit("row count mismatch")
    for covidx_line, bedtools_line, (x, y) in zip(
        covidx_rows, bedtools_rows, windows
    ):
        got = int(covidx_line.split("\t")[4])
        want = int(bedtools_line.split("\t")[4])
        if got != want:
            raise SystemExit(
                f"mismatch at {args.contig}:{x}-{y}: "
                f"covidx={got} bedtools={want}"
            )
    print(f"ok: {len(windows)} windows on {args.contig}")


if __name__ == "__main__":
    main()
