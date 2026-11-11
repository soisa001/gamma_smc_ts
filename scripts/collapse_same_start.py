#!/usr/bin/env python3
"""Truvari 5.5.0 collapse with an explicit same-CHROM/POS partition boundary.

Only candidate grouping differs from collapse_main. Matching, representative
selection, genotype consolidation and output annotations use native functions.
Run inside run_popped_simple_sv.py's pinned uv environment.
"""
from collections import defaultdict
import importlib.metadata
import logging
import sys

import truvari
import truvari.collapse as native


def same_start_chunks(chunks):
    count = 0
    for chunk, _ in chunks:
        groups = defaultdict(lambda: {"base": [], "__filtered": []})
        for category in ("base", "__filtered"):
            for record in chunk[category]:
                groups[record.chrom, record.pos][category].append(record)
        for locus in sorted(groups):
            count += 1
            if count % 10000 == 0:
                logging.info("SAME_START_PROGRESS groups=%d", count)
            yield groups[locus], count
    logging.info("SAME_START_COMPLETE groups=%d", count)


def main(argv=None):
    if importlib.metadata.version("truvari") != "5.5.0":
        raise ValueError("This adapter requires Truvari 5.5.0")
    args = native.parse_args(sys.argv[1:] if argv is None else argv)
    truvari.setup_logging(args.debug, show_version=True)
    if native.check_params(args):
        raise SystemExit("Invalid native collapse parameters")
    params = truvari.VariantParams(args=args, short_circuit=True, skip_gt=True,
                                  sizefilt=args.sizemin, chunksize=args.refdist)
    params.keep, params.hap, params.gt, params.chain = args.keep, args.hap, args.gt, args.chain
    params.sorter = native.SORTS[args.keep]
    params.consolidate = not args.no_consolidate
    base = truvari.VariantFile(args.input, params=params)
    regions = truvari.build_region_tree(base, includebed=args.bed)
    truvari.merge_region_tree_overlaps(regions)
    chunks = same_start_chunks(truvari.chunker(params, ('base', base.fetch_regions(regions))))
    chunks = native.tree_size_chunker(params, chunks)
    chunks = native.tree_dist_chunker(params, chunks, args.fast_cluster)
    outputs = native.CollapseOutput(args, params)
    try:
        for chunk in chunks:
            for call in native.collapse_chunk(chunk, params):
                # Assert the requested boundary on every emitted cluster.
                if any((m.comp.chrom, m.comp.pos) != (call.entry.chrom, call.entry.pos) for m in call.matches):
                    raise ValueError("Native cluster crossed a same-start partition")
                outputs.write(call, args.median_info)
    finally:
        outputs.close()
        base.close()
    outputs.dump_log()
    logging.info("Finished same-start collapse; native genotype consolidation")


if __name__ == "__main__":
    main()
