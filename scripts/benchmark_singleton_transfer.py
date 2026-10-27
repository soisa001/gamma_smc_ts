#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Deterministic synthetic parser/transfer benchmark, not a genomic accuracy study."""
import argparse
import json
from pathlib import Path
import random
import subprocess
import sys
import time

import phase_singletons as ph


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--samples", type=int, default=12000)
    parser.add_argument("--contigs", type=int, default=4)
    parser.add_argument("--anchors", type=int, default=40)
    parser.add_argument("--singletons", type=int, default=100, help="Per contig")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--chunks-per-chrom", type=int, default=1)
    parser.add_argument("--overlap-bp", type=int, default=1000000)
    args = parser.parse_args()
    if args.samples < 2 or args.contigs < 1 or args.anchors < 2 or args.singletons < 1:
        parser.error("Need >=2 samples/anchors and >=1 contig/singleton")
    args.out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(ph.SEED)
    samples = [f"synthetic_{i:05d}" for i in range(args.samples)]
    def coordinate(index):
        return index * 1000000 // (args.anchors + args.singletons + 1)
    jobs = []
    for chrom in range(1, args.contigs + 1):
        contig = f"chr{chrom}"
        header = (f"##fileformat=VCFv4.2\n##contig=<ID={contig},length=1000000>\n"
                  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                  '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
                  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        scaffold, hiphase = args.out / f"{contig}.scaffold.vcf", args.out / f"{contig}.hiphase.vcf"
        flips = [rng.randrange(2) for _ in samples]
        with open(scaffold, "w", newline="\n") as sf, open(hiphase, "w", newline="\n") as hf:
            sf.write(header)
            hf.write(header)
            left_count = args.anchors // 2
            for i in range(left_count):
                prefix = f"{contig}\t{coordinate(i+1)}\ta{i}\tA\tG\t.\tPASS\t.\tGT:PS\t"
                sf.write(prefix + "\t".join(["0|1:1"] * args.samples) + "\n")
                hf.write(prefix + "\t".join(f"{f}|{1-f}:1" for f in flips) + "\n")
            for i in range(args.singletons):
                carrier = rng.randrange(args.samples)
                calls = ["0|0:."] * args.samples
                # Every fifth singleton deliberately lacks a physical phase set.
                calls[carrier] = "0/1:." if i % 5 == 0 else f"{flips[carrier]}|{1-flips[carrier]}:1"
                hf.write(f"{contig}\t{coordinate(left_count+i+1)}\tr{i}\tC\tT\t.\tPASS\t.\tGT:PS\t" + "\t".join(calls) + "\n")
            for i in range(args.anchors - left_count):
                pos = coordinate(left_count + args.singletons + i + 1)
                prefix = f"{contig}\t{pos}\ta{left_count+i}\tA\tG\t.\tPASS\t.\tGT:PS\t"
                sf.write(prefix + "\t".join(["0|1:1"] * args.samples) + "\n")
                hf.write(prefix + "\t".join(f"{f}|{1-f}:1" for f in flips) + "\n")
        jobs.append({"contig": contig, "scaffold": str(scaffold.resolve()), "hiphase": str(hiphase.resolve())})
    manifest = args.out / "manifest.json"
    ph.json_write(manifest, {"genome_build": "synthetic", "scaffold_phase_scope": "chromosome", "samples": samples, "jobs": jobs})
    start = time.perf_counter()
    subprocess.run([sys.executable, str(Path(ph.__file__).resolve()), "--manifest", str(manifest), "--out", str(args.out / "results"), "--workers", str(args.workers),
                    "--chunks-per-chrom", str(args.chunks_per_chrom), "--overlap-bp", str(args.overlap_bp)], check=True)
    elapsed = time.perf_counter() - start
    summary = json.loads((args.out / "results" / "summary.json").read_text())
    # Assert synthetic truth, not just successful process exit.
    for result in summary["results"]:
        expected_excluded = (args.singletons + 4) // 5
        assert result["final_counts"] == {"physical": args.singletons - expected_excluded, "excluded": expected_excluded}
        with ph.text_input(result["output"]) as handle:
            for row, _ in ph.VCF(handle, result["contig"]).records():
                if ph.info_value(row, "SP_METHOD") == "PHYSICAL":
                    calls = ph.calls(row)
                    assert calls[ph.singleton(calls)][0] == (0, 1)
    input_bytes = sum(Path(job[k]).stat().st_size for job in jobs for k in ("scaffold", "hiphase"))
    report = {"kind": "synthetic correctness and parser throughput ONLY", "seed": ph.SEED,
              "parameters": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
              "wall_seconds": elapsed, "input_bytes_uncompressed": input_bytes,
              "input_MB_per_second": input_bytes / 1e6 / elapsed,
              "all_expected_phases_correct": True, "cached_jobs": sum(r["cached"] for r in summary["results"]),
              "final_counts": [r["final_counts"] for r in summary["results"]],
              "cross_chunk_lookups": sum(c["phase"].get("cross_chunk_lookups", 0) for r in summary["results"] for c in r.get("chunks", [])),
              "note": "Short synthetic GT:PS records; not real compressed-file, biological-accuracy, or full-cohort runtime evidence."}
    ph.json_write(args.out / "benchmark.json", report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
