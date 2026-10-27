#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Create an ordered cohort manifest without printing sample IDs or genotypes."""
import argparse
from pathlib import Path

import phase_singletons as ph


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scaffold", required=True, help="Local path, optionally containing {contig}")
    parser.add_argument("--hiphase", required=True, help="Local joint-cohort path, optionally containing {contig}")
    parser.add_argument("--contigs", required=True, nargs="+")
    parser.add_argument("--genome-build", required=True)
    parser.add_argument("--assert-chromosome-phase", required=True, action="store_true", help="Confirm the panel haplotypes are globally oriented within each chromosome")
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    samples, jobs = None, []
    for contig in args.contigs:
        job = {"contig": contig}
        for name in ("scaffold", "hiphase"):
            path = Path(getattr(args, name).format(contig=contig)).resolve()
            if not path.is_file():
                parser.error(f"Missing {name}: {path}")
            # Read the header only (not an indexed region) to avoid SIGPIPE
            # from an early-closed bcftools data stream.
            if path.suffix == ".bcf":
                import io
                import subprocess
                text = subprocess.check_output([args.bcftools, "view", "--no-version", "-h", str(path)], text=True)
                vcf_samples = ph.VCF(io.StringIO(text), contig).samples
            else:
                with ph.text_input(path, args.bcftools) as handle:
                    vcf_samples = ph.VCF(handle, contig).samples
            if samples is None:
                samples = vcf_samples
            elif samples != vcf_samples:
                parser.error(f"Cohort sample identity/order differs in {path}")
            job[name] = str(path)
        jobs.append(job)
    manifest = {"genome_build": args.genome_build, "scaffold_phase_scope": "chromosome", "samples": samples, "jobs": jobs}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    ph.json_write(args.output, manifest)
    print(f"Wrote {len(jobs)} contigs and {len(samples)} ordered samples to {args.output.resolve()}")


if __name__ == "__main__":
    main()
