#!/usr/bin/env python3
"""Native Truvari 5.5.0 collapse within original-POS groups, read sequentially.

Working VCFs are ordered by PM_ORIG_POS, not their trimmed POS, and intentionally
have no index. No reference alignment, chain extension, or GT patching is used.
"""
import csv
import importlib.metadata
from itertools import groupby
import logging
import sys

import truvari
import truvari.collapse as native


def main():
    if importlib.metadata.version('truvari') != '5.5.0':
        raise ValueError('Requires Truvari 5.5.0')
    args = native.parse_args(sys.argv[1:])
    if args.chain or args.hap or args.reference or args.bed or args.typeignore:
        raise ValueError('This adapter requires no chain/hap/reference/bed/typeignore')
    truvari.setup_logging(args.debug, show_version=True)
    params = truvari.VariantParams(args=args, short_circuit=True, skip_gt=True,
                                  sizefilt=args.sizemin, chunksize=args.refdist)
    params.keep, params.hap, params.gt, params.chain = args.keep, args.hap, args.gt, False
    params.sorter, params.consolidate = native.SORTS[args.keep], not args.no_consolidate
    base = truvari.VariantFile(args.input, params=params)
    output = native.CollapseOutput(args, params)
    previous = None
    try:
        with open(args.output+'.mapping.tsv', 'w') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['working_id', 'representative', 'cluster', 'members'])
            for number, (key, records) in enumerate(groupby(base, lambda r:(r.chrom, int(r.info['PM_ORIG_POS']))), 1):
                if previous is not None and key <= previous:
                    raise ValueError('Working records must be grouped in original chromosome/POS order')
                previous = key
                entries = list(records)
                kinds = {r.var_type() for r in entries}
                if len(kinds) != 1:
                    raise ValueError('Mixed types in an original-POS partition')
                for call in native.collapse_chunk(({'base':entries, '__filtered':[]}, number), params):
                    members = [call.entry, *[m.comp for m in call.matches]]
                    if any((r.chrom, int(r.info['PM_ORIG_POS'])) != key for r in members):
                        raise ValueError('Cluster crossed original-POS boundary')
                    for record in members:
                        writer.writerow([record.id, call.entry.id, call.match_id, len(members)])
                    output.write(call, args.median_info)
                if number % 10000 == 0:
                    logging.info('ORIGINAL_POSITION_PROGRESS groups=%d', number)
    finally:
        output.close()
        base.close()
    output.dump_log()


if __name__ == '__main__':
    main()
