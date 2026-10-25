#!/usr/bin/env python3
"""Compare covidx covered counts with samtools depth.

The mapping is documented in docs/parity.md. This script deliberately uses
only samtools and the covidx-query CLI, so it also works where pysam is absent.

Without --region or --windows-bed, the script enumerates short windows across
a whole contig and is intended for small contigs. Use --region for a bounded
read-rich region on a large contig. The region syntax is 1-based inclusive,
matching covidx-query --region.
"""
import argparse
import bisect
import json
import random
import re
import subprocess
import tempfile
from pathlib import Path

REGION_RE = re.compile(r"^(.+):(\d+)-(\d+)$")
MAX_WINDOWS = 2_000_000


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


def parse_region(text, contig):
    match = REGION_RE.match(text)
    if not match:
        raise SystemExit("--region must be CHROM:START-END")
    if match.group(1) != contig:
        raise SystemExit("--region contig disagrees with --contig")
    start = int(match.group(2)) - 1
    end = int(match.group(3))
    if start < 0 or end <= start:
        raise SystemExit("--region has invalid 1-based coordinates")
    return start, end


def short_window_count(length):
    if length <= 0:
        return 0
    if length <= 50:
        return length * (length + 1) // 2
    return 50 * length - 1225


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


def windows_for_region(start, end):
    windows = []
    for x in range(start, end):
        for size in range(1, min(50, end - x) + 1):
            windows.append((x, x + size))
    return windows


def load_windows(path, contig):
    windows = []
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith(("#", "track", "browser")):
            continue
        fields = line.split("\t")
        if fields[0] == contig:
            windows.append((int(fields[1]), int(fields[2])))
    return windows


def merge_windows(windows):
    merged = []
    for x, y in sorted(windows):
        if merged and x <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], y))
        else:
            merged.append((x, y))
    return merged


def write_regions(path, contig, intervals):
    path.write_text("".join(f"{contig}\t{x}\t{y}\n" for x, y in intervals))


def check_index_policy(index, args, ambiguous):
    path = Path(index) / "manifest.json"
    if not path.exists():
        raise SystemExit(f"cannot find manifest at {path}")
    policy = json.loads(path.read_text())["policy"]
    expected = {
        "min_mapq": args.min_mapq,
        "exclude_flags": args.exclude_flags,
        "include_supplementary": args.include_supplementary,
        "count_deletions": not args.no_count_deletions,
        "ambiguous_read_bases": ambiguous,
    }
    for key, value in expected.items():
        if policy.get(key) != value:
            raise SystemExit(
                f"index policy {key}={policy.get(key)!r} does not match "
                f"requested {value!r}"
            )


def resolve_read_group(index, sample, explicit):
    if explicit:
        return [part.strip() for part in explicit.split(",") if part.strip()]
    if not sample:
        return []
    path = Path(index) / "samples.tsv"
    if not path.exists():
        raise SystemExit(f"cannot resolve --sample without {path}")
    groups = []
    found = False
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 2 or fields[0] != sample:
            continue
        found = True
        if fields[1] != "-":
            groups.append(fields[1])
    if not found:
        raise SystemExit(f"sample {sample!r} is not listed in {path}")
    return groups


def prepare_bam(bam, read_group, regions, tmp):
    if not read_group:
        return bam
    filtered = Path(tmp) / "filtered.bam"
    command = ["samtools", "view", "-b"]
    for group in read_group:
        command += ["-r", group]
    command += ["-L", str(regions), "-o", str(filtered), bam]
    run(command)
    run(["samtools", "index", str(filtered)])
    return str(filtered)


def depth_positive(bam, contig, regions, sparse, args):
    effective = (args.exclude_flags | 0x800) & ~(
        0x800 if args.include_supplementary else 0
    )
    default = 0x704
    remove = default & ~effective
    add = effective & ~default
    command = [
        "samtools",
        "depth",
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
    if sparse:
        command += ["-b", str(regions)]
    else:
        command += ["-aa", "-r", contig]
    command.append(bam)
    positive = set()
    for line in run(command).stdout.splitlines():
        fields = line.split("\t")
        if len(fields) >= 3 and int(fields[2]) > 0:
            positive.add(int(fields[1]) - 1)
    return sorted(positive)


def covered_total(positive, windows):
    total = 0
    for x, y in windows:
        total += bisect.bisect_left(positive, y)
        total -= bisect.bisect_left(positive, x)
    return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--contig", required=True)
    parser.add_argument("--index", required=True)
    parser.add_argument("--unit", default=None)
    parser.add_argument("--sample", default=None)
    parser.add_argument("--read-group", default=None)
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
    parser.add_argument(
        "--region",
        default=None,
        help="1-based inclusive CHROM:START-END for bounded validation",
    )
    parser.add_argument("--windows-bed", default=None)
    parser.add_argument("--random-windows", type=int, default=500)
    parser.add_argument(
        "--allow-empty",
        action="store_true",
        help="permit a run with zero covered positions",
    )
    args = parser.parse_args()

    if args.sample and args.read_group:
        raise SystemExit("use either --sample or --read-group, not both")
    if args.ambiguous_read_bases == "exclude":
        raise SystemExit(
            "samtools depth cannot validate ambiguous-base exclusion; "
            "use naive_reference.py"
        )
    check_index_policy(args.index, args, args.ambiguous_read_bases)
    read_group = resolve_read_group(args.index, args.sample, args.read_group)
    if args.exclude_flags > 0xFFFF:
        raise SystemExit("--exclude-flags does not fit in 16 bits")
    if args.region and args.windows_bed:
        raise SystemExit("use either --region or --windows-bed, not both")
    length = contig_length(args.bam, args.contig)

    if args.region:
        start, end = parse_region(args.region, args.contig)
        if end > length:
            raise SystemExit("--region exceeds the contig length")
        if short_window_count(end - start) > MAX_WINDOWS:
            raise SystemExit(
                "region generates too many windows; use --windows-bed"
            )
        windows = windows_for_region(start, end)
        sparse = True
    elif args.windows_bed:
        windows = load_windows(args.windows_bed, args.contig)
        sparse = True
    else:
        if short_window_count(length) + args.random_windows > MAX_WINDOWS:
            raise SystemExit(
                "whole-contig windows exceed the limit; "
                "use --region or --windows-bed"
            )
        windows = windows_for(length, args.random_windows)
        sparse = False
    if not windows:
        raise SystemExit("no windows selected")
    for x, y in windows:
        if x < 0 or y <= x or y > length:
            raise SystemExit(f"invalid window {args.contig}:{x}-{y}")

    intervals = merge_windows(windows) if sparse else [(0, length)]
    with tempfile.TemporaryDirectory(prefix="covidx-samtools-") as tmp:
        regions = Path(tmp) / "regions.bed"
        write_regions(regions, args.contig, intervals)
        bam = prepare_bam(args.bam, read_group, regions, tmp)
        positive = depth_positive(bam, args.contig, regions, sparse, args)
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

    if covered_total(positive, windows) == 0 and not args.allow_empty:
        raise SystemExit(
            "no covered positions in the selected windows; parity would "
            "be vacuous (use --allow-empty for an empty case)"
        )

    rows = [
        line
        for line in query_output.splitlines()
        if line and not line.startswith("#")
    ]
    if len(rows) != len(windows):
        raise SystemExit(
            f"covidx-query returned {len(rows)} rows for {len(windows)} "
            "windows"
        )
    for line, (x, y) in zip(rows, windows):
        fields = line.split("\t")
        got = int(fields[4])
        want = bisect.bisect_left(positive, y)
        want -= bisect.bisect_left(positive, x)
        if got != want:
            raise SystemExit(
                f"mismatch at {args.contig}:{x}-{y}: "
                f"covidx={got} samtools={want}"
            )
    print(
        f"ok: {len(windows)} windows on {args.contig} "
        f"({len(positive)} positive positions)"
    )


if __name__ == "__main__":
    main()
