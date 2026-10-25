#!/usr/bin/env python3
"""Measured stress matrix for covidx build and query behavior.

Each stage varies one family of input-BAM or tool parameters and records wall
time, peak RSS, index size, run counts, and query throughput. Outputs JSON and
TSV summaries under the work directory.

Examples:
  ./bench/stress_matrix.py --stage synthetic
  ./bench/stress_matrix.py --stage threads --stage policies --stage query
  ./bench/stress_matrix.py --stage all --work-dir /tmp/covidx-stress
"""

import argparse
import hashlib
import json
import os
import random
import shutil
import signal
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEFAULT_BUILD = REPO / "build-release"


def run_measured(cmd, stdout_path, stderr_path, timeout=3600, env_extra=None):
    """Run cmd and return exit code, wall seconds, and peak RSS (MiB)."""
    env = os.environ.copy()
    if env_extra:
        env.update(env_extra)
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    pid = os.fork()
    if pid == 0:
        try:
            signal.alarm(timeout)
            with open(stdout_path, "wb") as out, open(stderr_path, "wb") as err:
                os.dup2(out.fileno(), 1)
                os.dup2(err.fileno(), 2)
                os.execvpe(cmd[0], cmd, env)
        except BaseException:
            os._exit(127)
    _, status, usage = os.wait4(pid, 0)
    elapsed = time.monotonic() - start
    code = os.waitstatus_to_exitcode(status)
    return {
        "exit": code,
        "wall_s": round(elapsed, 4),
        "maxrss_mb": round(usage.ru_maxrss / 1024.0, 1),
    }


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tree_hashes(root):
    hashes = {}
    for path in sorted(Path(root).rglob("*")):
        if path.is_file():
            hashes[str(path.relative_to(root))] = sha256_file(path)
    return hashes


def tree_size(root):
    return sum(
        path.stat().st_size for path in Path(root).rglob("*") if path.is_file()
    )


def manifest_for(index_dir):
    return json.loads((Path(index_dir) / "manifest.json").read_text())


def load_units(index_dir):
    manifest = manifest_for(index_dir)
    return [unit["name"] for unit in manifest["units"]]


def choose_unit(index_dir, preferred=None):
    names = load_units(index_dir)
    if preferred and preferred in names:
        return preferred
    if "pooled" in names:
        return "pooled"
    return names[0]


def write_windows_bed(path, contigs, count=5000, seed=1234):
    """contigs: list of (name, length). Deterministic mixed-width windows."""
    rng = random.Random(seed)
    widths = [1, 2, 10, 100, 1000, 10000, 100000, 1000000]
    lines = []
    usable = [(name, length) for name, length in contigs if length > 0]
    for i in range(count):
        name, length = usable[i % len(usable)]
        width = widths[rng.randrange(len(widths))]
        width = min(width, length)
        start = rng.randrange(0, length - width + 1)
        lines.append(f"{name}\t{start}\t{start + width}\n")
    # Guarantee boundary coverage for every contig (start, last base, and a
    # window ending exactly at the contig end).
    for name, length in usable:
        lines.append(f"{name}\t0\t1\n")
        lines.append(f"{name}\t{length - 1}\t{length}\n")
        lines.append(f"{name}\t0\t{length}\n")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("".join(lines))
    return len(lines)


def write_scope_windows_bed(path, manifest, count=5000, seed=1234):
    """Windows that stay inside one declared scope interval."""
    rng = random.Random(seed)
    widths = [1, 2, 10, 100, 1000, 10000, 100000]
    intervals = [
        (item["contig"], item["start"], item["end"])
        for item in manifest["scope"]
        if item["end"] > item["start"]
    ]
    lines = []
    for i in range(count):
        name, scope_start, scope_end = intervals[i % len(intervals)]
        width = min(widths[rng.randrange(len(widths))], scope_end - scope_start)
        start = rng.randrange(scope_start, scope_end - width + 1)
        lines.append(f"{name}\t{start}\t{start + width}\n")
    for name, scope_start, scope_end in intervals:
        lines.append(f"{name}\t{scope_start}\t{scope_start + 1}\n")
        lines.append(f"{name}\t{scope_end - 1}\t{scope_end}\n")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("".join(lines))
    return len(lines)


def count_rows(path):
    with open(path, "rb") as handle:
        return sum(1 for line in handle if line.strip())


class Stress:
    def __init__(self, args):
        self.args = args
        self.work = Path(args.work_dir)
        self.work.mkdir(parents=True, exist_ok=True)
        self.build_bin = str(Path(args.build_dir) / "covidx-build")
        self.query_bin = str(Path(args.query_dir) / "covidx-query")
        self.generator = str(Path(args.build_dir) / "covidx_stress_bam")
        self.results = []
        self.determinism = []

    def result(self, stage, **fields):
        row = {"stage": stage, **fields}
        self.results.append(row)
        print(json.dumps(row, sort_keys=True), flush=True)

    def generate(self, case):
        bam = self.work / "bams" / f"{case['name']}.bam"
        if bam.exists() and not self.args.regenerate:
            return bam
        bam.parent.mkdir(parents=True, exist_ok=True)
        cmd = [
            self.generator,
            "--out",
            str(bam),
            "--contigs",
            str(case.get("contigs", 4)),
            "--contig-length",
            str(case.get("length", 10_000_000)),
            "--reads-per-contig",
            str(case.get("reads", 200_000)),
            "--read-length",
            str(case.get("read_length", 150)),
            "--cigar",
            case.get("cigar", "m"),
            "--mapq",
            case.get("mapq", "constant"),
            "--flags",
            case.get("flags", "clean"),
            "--samples",
            str(case.get("samples", 1)),
            "--unmapped-reads",
            str(case.get("unmapped", 0)),
            "--cluster-percent",
            str(case.get("cluster", 0)),
            "--unknown-rg-percent",
            str(case.get("unknown_rg", 0)),
            "--no-rg-percent",
            str(case.get("no_rg", 0)),
            "--n-base-percent",
            str(case.get("n_base_percent", 0)),
            "--index",
            case.get("index", "auto"),
            "--seed",
            str(case.get("seed", 20261024)),
        ]
        measured = run_measured(
            cmd,
            bam.with_suffix(".gen.out"),
            bam.with_suffix(".gen.err"),
            timeout=3600,
        )
        if measured["exit"] != 0:
            raise RuntimeError(
                f"generator failed for {case['name']}: "
                f"{bam.with_suffix('.gen.err').read_text()[-2000:]}"
            )
        return bam

    def build(self, bam, case, threads=8, extra=None, tag=None):
        name = case["name"] if tag is None else tag
        out = self.work / "indexes" / name
        if out.exists():
            shutil.rmtree(out)
        out.mkdir(parents=True, exist_ok=True)
        work = self.work / "work" / name
        if work.exists():
            shutil.rmtree(work)
        cmd = [
            self.build_bin,
            "--bam",
            str(bam),
            "--out",
            str(out),
            "--threads",
            str(threads),
            "--work-dir",
            str(work),
        ]
        cmd += case.get("build_args", [])
        if extra:
            cmd += extra
        measured = run_measured(
            cmd,
            out / "build.out",
            out / "build.err",
            timeout=7200,
            env_extra={"SOURCE_DATE_EPOCH": "0"},
        )
        if measured["exit"] != 0:
            raise RuntimeError(
                f"build failed for {name}: "
                f"{(out / 'build.err').read_text()[-2000:]}"
            )
        index = out / f"{Path(bam).stem}.covidx"
        manifest = manifest_for(index)
        measured.update(
            {
                "index_dir": str(index),
                "index_mb": round(tree_size(index) / (1024 * 1024), 2),
                "runs": manifest["index"]["total_runs"],
                "covered": sum(
                    unit["total_covered"] for unit in manifest["units"]
                ),
                "units": len(manifest["units"]),
                "features": sum(
                    unit["total_features"] for unit in manifest["units"]
                ),
                "accepted": manifest["source"]["records_accepted"],
                "examined": manifest["source"]["records_examined"],
                "partial_scope": manifest["index"]["partial_scope"],
            }
        )
        return index, manifest, measured

    def query(
        self,
        index,
        unit,
        bed,
        validation="structural",
        crc=False,
        mmap=True,
        tag="query",
        extra=None,
    ):
        cmd = [
            self.query_bin,
            "--index",
            str(index),
            "--unit",
            unit,
            "--bed",
            str(bed),
            "--no-header",
            "--validation",
            validation,
            "--mmap" if mmap else "--no-mmap",
            "--no-validate-crc" if not crc else "--validate-crc",
        ]
        if extra:
            cmd += extra
        out = self.work / "query" / f"{tag}.tsv"
        measured = run_measured(cmd, out, out.with_suffix(".err"), timeout=3600)
        rows = count_rows(out) if measured["exit"] == 0 else 0
        measured.update({"rows": rows, "query_dir": str(out)})
        return measured

    def load_only(
        self, index, tag="load", validation="structural", crc=False, mmap=True
    ):
        cmd = [
            self.query_bin,
            "--index",
            str(index),
            "--list-units",
            "--validation",
            validation,
            "--mmap" if mmap else "--no-mmap",
            "--no-validate-crc" if not crc else "--validate-crc",
        ]
        out = self.work / "query" / f"{tag}.txt"
        return run_measured(cmd, out, out.with_suffix(".err"), timeout=3600)

    # ------------------------------------------------------------------
    # Stages
    # ------------------------------------------------------------------

    def stage_synthetic(self):
        cases = [
            {
                "name": "base_800k",
                "contigs": 4,
                "length": 10_000_000,
                "reads": 200_000,
            },
            {
                "name": "contigs_1",
                "contigs": 1,
                "length": 10_000_000,
                "reads": 800_000,
            },
            {
                "name": "contigs_64",
                "contigs": 64,
                "length": 1_000_000,
                "reads": 12_500,
            },
            {
                "name": "contigs_256",
                "contigs": 256,
                "length": 1_000_000,
                "reads": 3_125,
            },
            {
                "name": "length_1mb_dense",
                "contigs": 4,
                "length": 1_000_000,
                "reads": 200_000,
            },
            {
                "name": "length_500mb",
                "contigs": 4,
                "length": 500_000_000,
                "reads": 200_000,
            },
            {
                "name": "length_max_2gb",
                "contigs": 1,
                "length": 2_147_483_647,
                "reads": 50_000,
                "index": "csi",
            },
            {
                "name": "reads_10k",
                "contigs": 4,
                "length": 10_000_000,
                "reads": 2_500,
            },
            {
                "name": "reads_4m",
                "contigs": 4,
                "length": 10_000_000,
                "reads": 1_000_000,
                "build_threads": 20,
            },
            {
                "name": "readlen_50",
                "contigs": 4,
                "length": 10_000_000,
                "reads": 200_000,
                "read_length": 50,
            },
            {
                "name": "readlen_1000",
                "contigs": 4,
                "length": 10_000_000,
                "reads": 200_000,
                "read_length": 1000,
            },
            {
                "name": "readlen_10000",
                "contigs": 4,
                "length": 10_000_000,
                "reads": 20_000,
                "read_length": 10000,
            },
            {"name": "cigar_deletions", "cigar": "deletions"},
            {"name": "cigar_skips", "cigar": "skips"},
            {"name": "cigar_indels", "cigar": "indels"},
            {"name": "cigar_softclip", "cigar": "softclip"},
            {"name": "cigar_hardclip", "cigar": "hardclip"},
            {"name": "cigar_mixed", "cigar": "mixed"},
            {"name": "mapq_varied", "mapq": "varied"},
            {"name": "mapq_unavailable", "mapq": "unavailable"},
            {"name": "flags_mixed", "flags": "mixed"},
            {
                "name": "samples_8",
                "samples": 8,
                "build_args": ["--layout", "both"],
            },
            {
                "name": "samples_64",
                "samples": 64,
                "reads": 62_500,
                "build_args": ["--layout", "per-sample"],
            },
            {
                "name": "samples_256",
                "contigs": 2,
                "reads": 50_000,
                "samples": 256,
                "build_args": ["--layout", "per-sample"],
            },
            {
                "name": "contigs_4096",
                "contigs": 4096,
                "length": 100_000,
                "reads": 1,
            },
            {
                "name": "contigs_8192",
                "contigs": 8192,
                "length": 100_000,
                "reads": 1,
            },
            {
                "name": "single_contig_many_runs",
                "contigs": 1,
                "length": 500_000_000,
                "reads": 800_000,
                "build_threads": 8,
            },
            {
                "name": "cigar_manyops",
                "contigs": 4,
                "reads": 10_000,
                "read_length": 1000,
                "cigar": "manyops",
            },
            {
                "name": "unknown_rg",
                "samples": 4,
                "unknown_rg": 10,
                "no_rg": 5,
                "build_args": ["--unknown-sample", "unknown"],
            },
            {
                "name": "no_rg",
                "samples": 0,
                "build_args": ["--unknown-sample", "unknown"],
            },
            {"name": "clustered_90", "cluster": 90},
            {"name": "unmapped_tail", "unmapped": 100_000},
            {"name": "csi_10mb", "index": "csi"},
        ]
        if self.args.quick:
            cases = cases[:6]
        if self.args.only:
            wanted = set(self.args.only)
            cases = [case for case in cases if case["name"] in wanted]
        for case in cases:
            bam = self.generate(case)
            threads = case.get("build_threads", self.args.threads)
            index, manifest, built = self.build(bam, case, threads=threads)
            contigs = [(c["name"], c["length"]) for c in manifest["contigs"]]
            bed = self.work / "query" / f"windows_{case['name']}.bed"
            rows = write_windows_bed(bed, contigs, count=self.args.windows)
            unit = choose_unit(index)
            query = self.query(
                index, unit, bed, tag=f"synthetic_{case['name']}"
            )
            load = self.load_only(index, tag=f"synthetic_{case['name']}_load")
            bam_mb = round(bam.stat().st_size / (1024 * 1024), 2)
            self.result(
                "synthetic",
                case=case["name"],
                bam_mb=bam_mb,
                contigs=case.get("contigs", 4),
                contig_length=case.get("length", 10_000_000),
                reads=case.get("contigs", 4) * case.get("reads", 200_000)
                + case.get("unmapped", 0),
                read_length=case.get("read_length", 150),
                cigar=case.get("cigar", "m"),
                mapq=case.get("mapq", "constant"),
                flags=case.get("flags", "clean"),
                samples=case.get("samples", 1),
                index_type=case.get("index", "auto"),
                threads=threads,
                build_s=built["wall_s"],
                build_rss_mb=built["maxrss_mb"],
                index_mb=built["index_mb"],
                units=built["units"],
                runs=built["runs"],
                features=built["features"],
                accepted=built["accepted"],
                examined=built["examined"],
                load_s=load["wall_s"],
                load_rss_mb=load["maxrss_mb"],
                query_s=query["wall_s"],
                query_rows=query["rows"],
                ns_per_row=round(query["wall_s"] * 1e9 / max(rows, 1), 1),
                query_rss_mb=query["maxrss_mb"],
                query_unit=unit,
            )
            # A second, single-thread build must be byte-identical.
            if case["name"] in {
                "base_800k",
                "cigar_mixed",
                "samples_8",
                "clustered_90",
                "length_max_2gb",
            }:
                index1, _, built1 = self.build(
                    bam, case, threads=1, tag=case["name"] + "_t1"
                )
                h1 = tree_hashes(index)
                h8 = tree_hashes(index1)
                self.determinism.append(
                    {
                        "stage": "synthetic",
                        "case": case["name"],
                        "identical": h1 == h8,
                        "files": len(h1),
                        "t1_s": built1["wall_s"],
                        "t8_s": built["wall_s"],
                    }
                )

    def stage_threads(self):
        bam = REPO / "data-raw" / "multi-sample.bam"
        if not bam.exists():
            print(
                "multi-sample.bam missing; skipping thread scaling",
                file=sys.stderr,
            )
            return
        case = {"name": "threads_multisample"}
        levels = [1, 2, 4, 8, 16, 30]
        if self.args.quick:
            levels = [1, 30]
        hashes = {}
        for level in levels:
            index, manifest, built = self.build(
                bam,
                case,
                threads=level,
                extra=["--layout", "both", "--overwrite"],
                tag=f"threads_t{level}",
            )
            hashes[level] = tree_hashes(index)
            bed = self.work / "query" / "threads_windows.bed"
            if not bed.exists():
                contigs = [
                    (c["name"], c["length"]) for c in manifest["contigs"]
                ]
                write_windows_bed(bed, contigs, count=self.args.windows)
            query = self.query(
                index, choose_unit(index), bed, tag=f"threads_t{level}"
            )
            self.result(
                "threads",
                threads=level,
                bam="multi-sample.bam",
                bam_mb=round(bam.stat().st_size / (1024 * 1024), 1),
                build_s=built["wall_s"],
                build_rss_mb=built["maxrss_mb"],
                index_mb=built["index_mb"],
                runs=built["runs"],
                features=built["features"],
                accepted=built["accepted"],
                units=built["units"],
                query_s=query["wall_s"],
                query_rows=query["rows"],
            )
        base = levels[0]
        for level in levels[1:]:
            self.determinism.append(
                {
                    "stage": "threads",
                    "case": f"threads t1 vs t{level}",
                    "identical": hashes[base] == hashes[level],
                    "files": len(hashes[level]),
                }
            )

    def stage_policies(self):
        bam = REPO / "data-raw" / "cs_t2t_se_cereale.bam"
        if not bam.exists():
            print(
                "cs_t2t_se_cereale.bam missing; skipping policies",
                file=sys.stderr,
            )
            return
        variants = [
            ("pooled_default", ["--layout", "pooled"]),
            ("per_sample", ["--layout", "per-sample"]),
            ("both", ["--layout", "both"]),
            ("min_mapq_20", ["--layout", "pooled", "--min-mapq", "20"]),
            ("min_mapq_30", ["--layout", "pooled", "--min-mapq", "30"]),
            ("min_mapq_60", ["--layout", "pooled", "--min-mapq", "60"]),
            ("min_mapq_255", ["--layout", "pooled", "--min-mapq", "255"]),
            ("flags_none", ["--layout", "pooled", "--exclude-flags", "0"]),
            (
                "flags_plus_dups",
                ["--layout", "pooled", "--exclude-flags", "0x304"],
            ),
            (
                "include_supplementary",
                ["--layout", "pooled", "--include-supplementary"],
            ),
            ("no_deletions", ["--layout", "pooled", "--no-count-deletions"]),
            (
                "ambiguous_exclude",
                ["--layout", "pooled", "--ambiguous-read-bases", "exclude"],
            ),
            ("checksum_source", ["--layout", "pooled", "--checksum-source"]),
            ("contigs_one", ["--layout", "pooled", "--contigs", "CP155614.1"]),
            (
                "scope_random",
                [
                    "--layout",
                    "pooled",
                    "--regions",
                    str(self.work / "bams" / "scope_cs.bed"),
                ],
            ),
        ]
        # Deterministic 80-interval scope inside the first 5 Mb of CP155614.1.
        scope_bed = self.work / "bams" / "scope_cs.bed"
        scope_bed.parent.mkdir(parents=True, exist_ok=True)
        rng = random.Random(77)
        intervals = []
        pos = 0
        while pos < 5_000_000 and len(intervals) < 80:
            pos += rng.randrange(0, 30_000)
            length = rng.randrange(1, 4_000)
            if pos + length > 5_000_000:
                break
            intervals.append((pos, pos + length))
            pos += length + rng.randrange(0, 20_000)
        scope_bed.write_text(
            "".join(f"CP155614.1\t{a}\t{b}\n" for a, b in intervals)
        )
        for name, extra in variants:
            case = {"name": f"policy_{name}"}
            index, manifest, built = self.build(
                bam,
                case,
                threads=self.args.threads,
                extra=extra + ["--overwrite"],
                tag=f"policy_{name}",
            )
            bed = self.work / "query" / f"policy_{name}.bed"
            if manifest["index"]["partial_scope"]:
                write_scope_windows_bed(
                    bed, manifest, count=self.args.windows, seed=99
                )
            else:
                contigs = [
                    (c["name"], c["length"]) for c in manifest["contigs"]
                ]
                write_windows_bed(
                    bed, contigs, count=self.args.windows, seed=99
                )
            unit = choose_unit(index)
            query = self.query(index, unit, bed, tag=f"policy_{name}")
            self.result(
                "policies",
                variant=name,
                bam="cs_t2t_se_cereale.bam",
                build_s=built["wall_s"],
                build_rss_mb=built["maxrss_mb"],
                index_mb=built["index_mb"],
                runs=built["runs"],
                features=built["features"],
                covered=built["covered"],
                accepted=built["accepted"],
                examined=built["examined"],
                units=built["units"],
                partial_scope=built["partial_scope"],
                query_s=query["wall_s"],
                query_rows=query["rows"],
            )

    def stage_query(self):
        multi_bam = REPO / "data-raw" / "multi-sample.bam"
        cs_bam = REPO / "data-raw" / "cs_t2t_se_cereale.bam"
        if not multi_bam.exists() or not cs_bam.exists():
            print("real BAMs missing; skipping query stage", file=sys.stderr)
            return
        multi_case = {"name": "query_multisample"}
        multi_index, manifest, _ = self.build(
            multi_bam,
            multi_case,
            threads=self.args.threads,
            extra=["--layout", "both", "--overwrite"],
            tag="query_multisample",
        )
        cs_case = {"name": "query_cs"}
        cs_index, cs_manifest, _ = self.build(
            cs_bam,
            cs_case,
            threads=self.args.threads,
            extra=["--layout", "pooled", "--overwrite"],
            tag="query_cs",
        )

        contigs = [(c["name"], c["length"]) for c in manifest["contigs"]]
        load_bed = self.work / "query" / "query_load.bed"
        write_windows_bed(load_bed, contigs, count=self.args.windows, seed=555)
        unit = choose_unit(multi_index)
        variants = [
            (
                "mmap_crc_structural",
                dict(mmap=True, crc=True, validation="structural"),
            ),
            (
                "mmap_nocrc_structural",
                dict(mmap=True, crc=False, validation="structural"),
            ),
            (
                "heap_nocrc_structural",
                dict(mmap=False, crc=False, validation="structural"),
            ),
            ("mmap_nocrc_none", dict(mmap=True, crc=False, validation="none")),
            ("mmap_nocrc_full", dict(mmap=True, crc=False, validation="full")),
        ]
        for name, options in variants:
            load = self.load_only(multi_index, tag=f"qload_{name}", **options)
            query = self.query(
                multi_index, unit, load_bed, tag=f"q_{name}", **options
            )
            self.result(
                "query_load",
                variant=name,
                index="multi-sample both",
                index_mb=round(tree_size(multi_index) / (1024 * 1024), 2),
                load_s=load["wall_s"],
                load_rss_mb=load["maxrss_mb"],
                query_s=query["wall_s"],
                rows=query["rows"],
                ns_per_row=round(
                    query["wall_s"] * 1e9 / max(query["rows"], 1), 1
                ),
                query_rss_mb=query["maxrss_mb"],
            )

        shapes = [
            ("bed_1k", 1_000),
            ("bed_10k", 10_000),
            ("bed_100k", 100_000),
            ("bed_500k", 500_000),
        ]
        for name, count in shapes:
            bed = self.work / "query" / f"shape_{name}.bed"
            write_windows_bed(bed, contigs, count=count, seed=1000 + count)
            query = self.query(multi_index, unit, bed, tag=f"shape_{name}")
            self.result(
                "query_shape",
                variant=name,
                index="multi-sample both",
                rows=query["rows"],
                query_s=query["wall_s"],
                ns_per_row=round(
                    query["wall_s"] * 1e9 / max(query["rows"], 1), 1
                ),
                query_rss_mb=query["maxrss_mb"],
            )

        # Tiled queries: row count and throughput for several window sizes.
        # The 1 kb tiling is restricted to one contig; tiling all 21 Gb at
        # 1 kb would emit 21M rows and measure I/O rather than query cost.
        tiles = [
            ("1kb_one_contig", 1_000, ["--contig", "CP155614.1"]),
            ("1mb_all", 1_000_000, []),
            ("10mb_all", 10_000_000, []),
            ("100mb_all", 100_000_000, []),
        ]
        for name, window, extra in tiles:
            out = self.work / "query" / f"tile_{name}.tsv"
            cmd = [
                self.query_bin,
                "--index",
                str(multi_index),
                "--unit",
                unit,
                "--window",
                str(window),
                "--no-header",
                "--no-validate-crc",
            ] + extra
            measured = run_measured(
                cmd, out, out.with_suffix(".err"), timeout=3600
            )
            rows = count_rows(out) if measured["exit"] == 0 else 0
            self.result(
                "query_tile",
                variant=name,
                window=window,
                index="multi-sample both",
                rows=rows,
                query_s=measured["wall_s"],
                ns_per_row=round(measured["wall_s"] * 1e9 / max(rows, 1), 1),
                query_rss_mb=measured["maxrss_mb"],
            )

        # Compare index size / run density across the two real BAMs.
        for label, index in (
            ("cs pooled", cs_index),
            ("multi both", multi_index),
        ):
            bed = self.work / "query" / f"density_{label.replace(' ', '_')}.bed"
            manifest_here = manifest_for(index)
            write_windows_bed(
                bed,
                [(c["name"], c["length"]) for c in manifest_here["contigs"]],
                count=self.args.windows,
                seed=777,
            )
            query = self.query(
                index, choose_unit(index), bed, tag=f"density_{label}"
            )
            load = self.load_only(index, tag=f"density_load_{label}")
            self.result(
                "query_density",
                index=label,
                index_mb=round(tree_size(index) / (1024 * 1024), 2),
                runs=manifest_here["index"]["total_runs"],
                features=sum(
                    unit["total_features"]
                    for unit in manifest_here["units"]
                ),
                load_s=load["wall_s"],
                load_rss_mb=load["maxrss_mb"],
                query_s=query["wall_s"],
                rows=query["rows"],
                ns_per_row=round(
                    query["wall_s"] * 1e9 / max(query["rows"], 1), 1
                ),
            )

        # Feature queries: pooled summation and the large-input batch path.
        # Every query emits feature counts, so there is no opt-in variant.
        feature_manifest = manifest_for(multi_index)
        sample_units = [
            item["name"]
            for item in feature_manifest["units"]
            if item["kind"] == "sample"
        ]
        feature_bed = self.work / "query" / "feature_windows.bed"
        write_windows_bed(
            feature_bed, contigs, count=min(self.args.windows, 2000), seed=888
        )
        rows_by_unit = {}
        for name in sample_units + ["pooled"]:
            query = self.query(
                multi_index,
                name,
                feature_bed,
                tag=f"feature_sum_{name}",
            )
            rows_by_unit[name] = [
                line.split("\t")
                for line in Path(query["query_dir"]).read_text().splitlines()
            ]
        mismatches = 0
        for index, row in enumerate(rows_by_unit["pooled"]):
            member_sum = sum(
                int(rows_by_unit[name][index][4]) for name in sample_units
            )
            if int(row[4]) != member_sum:
                mismatches += 1
        self.result(
            "query_features",
            variant="pooled_sum",
            index="multi-sample both",
            units=len(sample_units),
            rows=len(rows_by_unit["pooled"]),
            mismatches=mismatches,
        )
        full_load = self.load_only(
            multi_index, tag="feature_full_load", validation="full"
        )
        self.result(
            "query_features",
            variant="full_validation",
            index="multi-sample both",
            load_s=full_load["wall_s"],
            load_rss_mb=full_load["maxrss_mb"],
        )
        large_bed = self.work / "query" / "shape_bed_500k.bed"
        if large_bed.exists():
            large = self.query(
                multi_index,
                unit,
                large_bed,
                tag="feature_500k",
            )
            self.result(
                "query_features",
                variant="large_500k",
                index="multi-sample both",
                rows=large["rows"],
                query_s=large["wall_s"],
                ns_per_row=round(
                    large["wall_s"] * 1e9 / max(large["rows"], 1), 1
                ),
                query_rss_mb=large["maxrss_mb"],
            )

    def stage_reference(self):
        """Independent samtools-based reference for CIGAR/policy variants."""
        reference_script = REPO / "data-raw" / "naive_reference.py"
        if not reference_script.exists():
            print(
                "naive_reference.py missing; skipping reference stage",
                file=sys.stderr,
            )
            return
        base = {
            "contigs": 2,
            "length": 300_000,
            "reads": 10_000,
            "read_length": 100,
            "samples": 1,
            "mapq": "constant",
            "flags": "clean",
        }
        cases = []
        for cigar in (
            "m",
            "deletions",
            "skips",
            "indels",
            "softclip",
            "hardclip",
            "eqx",
            "mixed",
        ):
            cases.append(
                ({**base, "name": f"ref_{cigar}", "cigar": cigar}, [], [])
            )
        cases.append(
            (
                {**base, "name": "ref_del_gap", "cigar": "deletions"},
                ["--no-count-deletions"],
                ["--no-count-deletions"],
            )
        )
        cases.append(
            (
                {
                    **base,
                    "name": "ref_n_cover",
                    "cigar": "m",
                    "n_base_percent": 10,
                },
                [],
                [],
            )
        )
        cases.append(
            (
                {
                    **base,
                    "name": "ref_n_exclude",
                    "cigar": "m",
                    "n_base_percent": 10,
                },
                ["--ambiguous-read-bases", "exclude"],
                ["--ambiguous-read-bases", "exclude"],
            )
        )
        cases.append(
            (
                {
                    **base,
                    "name": "ref_flags_mapq",
                    "cigar": "mixed",
                    "flags": "mixed",
                    "mapq": "varied",
                },
                ["--min-mapq", "20"],
                ["--min-mapq", "20"],
            )
        )
        cases.append(
            (
                {
                    **base,
                    "name": "ref_include_supp",
                    "cigar": "mixed",
                    "flags": "mixed",
                },
                ["--include-supplementary"],
                ["--include-supplementary"],
            )
        )
        for case, build_args, reference_args in cases:
            bam = self.generate(case)
            index, manifest, built = self.build(
                bam,
                case,
                threads=4,
                extra=["--layout", "pooled"] + build_args,
                tag=case["name"],
            )
            contigs = [(c["name"], c["length"]) for c in manifest["contigs"]]
            bed = self.work / "reference" / f"{case['name']}.bed"
            write_windows_bed(bed, contigs, count=500, seed=4242)
            cmd = [
                sys.executable,
                str(reference_script),
                "--bam",
                str(bam),
                "--contig",
                "chr1",
                "--index",
                str(index),
                "--unit",
                "pooled",
                "--covidx-query",
                self.query_bin,
                "--windows-bed",
                str(bed),
            ] + reference_args
            out = self.work / "reference" / f"{case['name']}.out"
            measured = run_measured(
                cmd, out, out.with_suffix(".err"), timeout=3600
            )
            self.result(
                "reference",
                case=case["name"],
                exit=measured["exit"],
                ok=measured["exit"] == 0,
                build_s=built["wall_s"],
                runs=built["runs"],
                features=built["features"],
                bam_mb=round(bam.stat().st_size / (1024 * 1024), 2),
                detail=(
                    out.read_text().strip()[-160:]
                    if measured["exit"] == 0
                    else out.with_suffix(".err").read_text().strip()[-240:]
                ),
            )

    def stage_correctness(self):
        """Cross-check synthetic coverage against the generator-independent
        full-validation load and check query rows stay in range."""
        index_root = self.work / "indexes"
        if not index_root.exists():
            return
        for index in sorted(index_root.glob("*/*.covidx")):
            if "t1" in index.parent.name or "policy_" in index.parent.name:
                continue
            manifest = manifest_for(index)
            unit = choose_unit(index)
            out = self.work / "query" / f"correct_{index.parent.name}.tsv"
            cmd = [
                self.query_bin,
                "--index",
                str(index),
                "--unit",
                unit,
                "--window",
                "1000000",
                "--no-header",
                "--validation",
                "full",
                "--no-validate-crc",
            ]
            measured = run_measured(
                cmd, out, out.with_suffix(".err"), timeout=3600
            )
            ok = measured["exit"] == 0
            bad = 0
            bad_features = 0
            rows = 0
            if ok:
                for line in out.read_text().splitlines():
                    fields = line.split("\t")
                    if len(fields) < 8:
                        bad += 1
                        bad_features += 1
                        continue
                    rows += 1
                    start, end = int(fields[1]), int(fields[2])
                    features, covered = int(fields[4]), int(fields[5])
                    if end <= start or covered > end - start or covered < 0:
                        bad += 1
                    if features < 0:
                        bad_features += 1
            self.result(
                "correctness",
                index=index.parent.name,
                exit=measured["exit"],
                rows=rows,
                bad_rows=bad,
                query_s=measured["wall_s"],
                feature_exit=measured["exit"],
                feature_rows=rows,
                bad_feature_rows=bad_features,
                feature_query_s=measured["wall_s"],
                full_validation_ok=ok and bad_features == 0,
            )

    @staticmethod
    def row_key(row):
        stage = row["stage"]
        if stage in ("synthetic", "reference"):
            identity = row.get("case", "")
        elif stage == "correctness":
            identity = row.get("index", "")
        elif "variant" in row:
            identity = row.get("variant", "")
        elif "threads" in row:
            identity = f"threads={row['threads']}"
        elif "index" in row:
            identity = row.get("index", "")
        else:
            identity = ""
        return f"{stage}\x1f{identity}"

    def write_outputs(self):
        # Preserve rows from earlier stage runs in the same work directory.
        # A rerun replaces the rows with the same stage and identity.
        combined = {}
        det_combined = {}
        results_path = self.work / "stress-results.json"
        if results_path.exists():
            try:
                previous = json.loads(results_path.read_text())
                for row in previous.get("results", []):
                    combined[self.row_key(row)] = row
                for row in previous.get("determinism", []):
                    row = dict(row)
                    if not row.get("stage"):
                        row["stage"] = (
                            "threads"
                            if str(row.get("case", "")).startswith("threads ")
                            else "synthetic"
                        )
                    det_combined[
                        f"{row.get('stage', '')}\x1f{row.get('case', '')}"
                    ] = row
            except (json.JSONDecodeError, KeyError):
                pass
        for row in self.results:
            combined[self.row_key(row)] = row
        for row in self.determinism:
            det_combined[f"{row.get('stage', '')}\x1f{row.get('case', '')}"] = (
                row
            )
        results = list(combined.values())
        determinism = list(det_combined.values())
        results_path.write_text(
            json.dumps(
                {
                    "results": results,
                    "determinism": determinism,
                },
                indent=2,
            )
            + "\n"
        )
        for stage in sorted({row["stage"] for row in results}):
            rows = [row for row in results if row["stage"] == stage]
            keys = sorted({key for row in rows for key in row})
            tsv = self.work / f"{stage}.tsv"
            with open(tsv, "w") as handle:
                handle.write("\t".join(keys) + "\n")
                for row in rows:
                    handle.write(
                        "\t".join(str(row.get(key, "")) for key in keys) + "\n"
                    )
        if determinism:
            det = self.work / "determinism.tsv"
            with open(det, "w") as handle:
                handle.write("stage\tcase\tidentical\tfiles\tt1_s\tt8_s\n")
                for row in determinism:
                    values = []
                    for key in (
                        "stage",
                        "case",
                        "identical",
                        "files",
                        "t1_s",
                        "t8_s",
                    ):
                        value = row.get(key)
                        values.append(
                            "-" if value is None or value == "" else str(value)
                        )
                    handle.write("\t".join(values) + "\n")
        print(f"\nwrote {results_path} ({len(results)} rows)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default="/tmp/covidx-stress")
    parser.add_argument("--build-dir", default=str(DEFAULT_BUILD))
    parser.add_argument(
        "--query-dir",
        default=None,
        help="directory holding covidx-query (default: --build-dir)",
    )
    parser.add_argument(
        "--stage",
        action="append",
        choices=[
            "synthetic",
            "threads",
            "policies",
            "query",
            "reference",
            "correctness",
            "all",
        ],
        default=None,
    )
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--windows", type=int, default=5000)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--regenerate", action="store_true")
    parser.add_argument(
        "--only",
        action="append",
        default=None,
        help="run only these synthetic case names (repeatable)",
    )
    args = parser.parse_args()
    if args.query_dir is None:
        args.query_dir = args.build_dir
    stages = args.stage or ["all"]
    if "all" in stages:
        stages = [
            "synthetic",
            "threads",
            "policies",
            "query",
            "reference",
            "correctness",
        ]
    stress = Stress(args)
    for stage in stages:
        print(f"=== stage: {stage} ===", flush=True)
        getattr(stress, f"stage_{stage}")()
    stress.write_outputs()


if __name__ == "__main__":
    main()
