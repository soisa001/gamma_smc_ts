#!/usr/bin/env python3
"""Truvari 5.5.0 clustering with phase-preserving --gt all consolidation.

Native matching, representative choice, -o and -c writers are unchanged.
Only gt_aware_consolidate is replaced: a sample's unique ALT-bearing member
supplies its GT, phasing flag and FORMAT fields. No haplotype is inferred.
"""
from __future__ import annotations

from collections import Counter
import logging
import sys

AUDIT = Counter()
FREQUENCY_HEADERS = (
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="ALT count after phased consolidation">',
    '##INFO=<ID=AN,Number=1,Type=Integer,Description="Called allele count after phased consolidation">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="AC divided by AN after phased consolidation">',
)


def ensure_frequency_headers(header):
    for key, line in zip(("AC", "AN", "AF"), FREQUENCY_HEADERS):
        if key not in header.info:
            header.add_line(line)


def copy_call(destination, source):
    """Copy the donor's complete FORMAT call, clearing unrelated old fields."""
    phase = source.phased
    donor = dict(source.items())
    for key in list(destination.keys()):
        if key not in donor:
            old = destination[key]
            destination[key] = tuple(None for _ in old) if isinstance(old, tuple) else None
    for key, value in donor.items():
        destination[key] = value
    destination.phased = phase
    if destination["GT"] != donor["GT"] or destination.phased != phase:
        raise ValueError("Donor GT/phase changed during consolidation")
    for key in ("PS", "PID", "PGT", "HP"):
        if key in donor and destination[key] != donor[key]:
            raise ValueError(f"Donor phase field {key} changed during consolidation")


def consolidate(entry, others):
    if not others:
        return entry, 0
    copied = ac = an = 0
    for sample in entry.samples:
        destination = entry.samples[sample]
        present = 1 in destination["GT"]
        donor = None
        for match in others:
            call = match.comp.samples[sample]
            if 1 in call["GT"]:
                if present or donor is not None:
                    raise ValueError("Multiple ALT-bearing members in one sample despite --gt all")
                donor = call
        if donor is not None:
            copy_call(destination, donor)
            copied += 1
            AUDIT["copied_phased_calls" if donor.phased else "copied_unphased_calls"] += 1
        # No ALT donor: preserve the representative's reference/missing call.
        gt = destination["GT"]
        ac += gt.count(1)
        an += sum(allele is not None for allele in gt)
    ensure_frequency_headers(entry.header)
    entry.info["AC"] = (ac,)
    entry.info["AN"] = an
    entry.info["AF"] = (ac / an if an else None,)
    AUDIT["merged_clusters"] += 1
    return entry, copied


def main(argv=None):
    import importlib.metadata
    import truvari.collapse as collapse
    argv = list(sys.argv[1:] if argv is None else argv)
    if importlib.metadata.version("truvari") != "5.5.0":
        raise SystemExit("This phase-preserving wrapper requires Truvari 5.5.0")
    parsed = collapse.parse_args(argv)
    if parsed.gt != "all" or parsed.hap or parsed.no_consolidate:
        raise SystemExit("Phase-preserving consolidation requires --gt all, without --hap or --no-consolidate")
    original_edit_header = collapse.edit_header
    def edit_header(*args, **kwargs):
        header = original_edit_header(*args, **kwargs)
        ensure_frequency_headers(header)
        header.add_line('##phased_consolidation="Unique ALT donor GT, phase flag and FORMAT copied; gt=all"')
        return header
    AUDIT.clear()
    collapse.edit_header = edit_header
    collapse.gt_aware_consolidate = consolidate
    collapse.collapse_main(argv)
    logging.info("PHASE_PRESERVATION_AUDIT %s", dict(AUDIT))


if __name__ == "__main__":
    main()
