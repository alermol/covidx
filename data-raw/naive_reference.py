#!/usr/bin/env python3
"""Naive per-base reference for covidx, including ambiguous-base exclusion.

Uses `samtools view` and plain Python; no pysam. This is the reference for
--ambiguous-read-bases exclude, where samtools depth is not a valid comparator.
"""
import argparse
import random
import re
import subprocess
from pathlib import Path

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=XB])")


def run(cmd):
    return subprocess.run(cmd, check=True, text=True, capture_output=True)


def contig_length(bam, contig):
    for line in run(["samtools", "view", "-H", bam]).stdout.splitlines():
        if line.startswith("@SQ"):
            fields = dict(part.split(":", 1) for part in line.split("\t")[1:])
            if fields.get("SN") == contig:
                return int(fields["LN"])
    raise SystemExit(f"unknown contig {contig!r}")


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
    parser.add_argument(
        "--ambiguous-read-bases", choices=["cover", "exclude"], default="cover"
    )
    parser.add_argument("--windows-bed", default=None)
    parser.add_argument("--random-windows", type=int, default=500)
    args = parser.parse_args()

    length = contig_length(args.bam, args.contig)
    rng = random.Random(31)
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
    covered = [False] * length
    records = run(
        ["samtools", "view", args.bam, args.contig]
    ).stdout.splitlines()
    for record in records:
        fields = record.split("\t")
        flag = int(fields[1])
        if flag & 4 or flag & effective:
            continue
        mapq = int(fields[4])
        if mapq < args.min_mapq:
            continue
        pos = int(fields[3]) - 1
        cigar = fields[5]
        seq = fields[9]
        if cigar == "*":
            raise SystemExit(f"mapped read {fields[0]} has CIGAR *")
        if args.ambiguous_read_bases == "exclude" and seq == "*":
            raise SystemExit(
                f"read {fields[0]} has SEQ=* but exclusion was requested"
            )
        ref = pos
        query = 0
        for count_text, op in CIGAR_RE.findall(cigar):
            count = int(count_text)
            if count == 0:
                raise SystemExit("zero-length CIGAR operation")
            if op in ("M", "=", "X"):
                for k in range(count):
                    include = True
                    if args.ambiguous_read_bases == "exclude":
                        include = seq[query + k].upper() in "ACGT"
                    if include:
                        covered[ref + k] = True
                ref += count
                query += count
            elif op == "D":
                if not args.no_count_deletions:
                    for k in range(count):
                        covered[ref + k] = True
                ref += count
            elif op == "N":
                ref += count
            elif op in ("I", "S"):
                query += count
            elif op in ("H", "P"):
                pass
            elif op == "B":
                raise SystemExit("unsupported CIGAR B")
            else:
                raise SystemExit(f"unsupported CIGAR operation {op}")

    prefix = [0] * (length + 1)
    for i, value in enumerate(covered):
        prefix[i + 1] = prefix[i] + (1 if value else 0)

    with __import__("tempfile").TemporaryDirectory(
        prefix="covidx-naive-"
    ) as tmp:
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
        output = run(query).stdout

    rows = [
        line
        for line in output.splitlines()
        if line and not line.startswith("#")
    ]
    if len(rows) != len(windows):
        raise SystemExit("row count mismatch")
    for line, (x, y) in zip(rows, windows):
        got = int(line.split("\t")[4])
        want = prefix[y] - prefix[x]
        if got != want:
            raise SystemExit(
                f"mismatch at {args.contig}:{x}-{y}: covidx={got} naive={want}"
            )
    print(f"ok: {len(windows)} windows on {args.contig}")


if __name__ == "__main__":
    main()
