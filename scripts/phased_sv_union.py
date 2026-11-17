"""Explicit haplotype union for an already chosen biallelic allele cluster.

No clustering, realignment, read evidence aggregation, or genotype imputation.
Mixed/unphased ambiguous calls and contradictory phase sets are flagged rather
than assigning ALT copies to guessed haplotypes.
"""


def snapshot(sample):
    return dict(gt=tuple(sample.get('GT') or ()),phased=sample.phased,ps=sample.get('PS'))


def union_calls(calls):
    if not calls or any(len(c['gt'])!=2 or any(x not in (0,1,None) for x in c['gt']) for c in calls):
        return dict(gt=(None,None),phased=False,ps=None,status='unsupported_gt')
    alt_calls = [c for c in calls if 1 in c['gt']]
    phase_sets = {c['ps'] for c in alt_calls if c['phased'] and c['ps'] is not None}
    if len(phase_sets)>1:
        return dict(gt=(None,None),phased=False,ps=None,status='phase_set_conflict')
    uncertain = [c for c in calls if not c['phased'] and c['gt'][0]!=c['gt'][1]]
    if uncertain:
        # Preserve an unphased donor's dosage when every other record is fully
        # reference. With multiple possible carriers, dosage itself can be ambiguous.
        if len(uncertain)==1 and all(c is uncertain[0] or c['gt']==(0,0) for c in calls):
            return dict(gt=uncertain[0]['gt'],phased=False,ps=None,status='unphased_single_donor')
        if any(c['gt']==(1,1) for c in calls):
            return dict(gt=(1,1),phased=True,ps=None,status='homozygous_alt_covers_both')
        return dict(gt=(None,None),phased=False,ps=None,status='unphased_union_ambiguous')
    result = []
    for hap in (0,1):
        values = [c['gt'][hap] for c in calls]
        result.append(1 if 1 in values else None if None in values else 0)
    ps = next(iter(phase_sets)) if len(phase_sets)==1 and all(c['ps'] is not None for c in alt_calls if c['phased']) else None
    return dict(gt=tuple(result),phased=True,ps=ps,
                status='partial_missing_union' if None in result else 'phased_union')


def apply_union(record, originals):
    """Set only GT/phase/PS on a *new trial record*, returning sample statuses."""
    statuses = {}
    for sample in record.samples:
        result = union_calls([snapshot(r.samples[sample]) for r in originals])
        record.samples[sample]['GT'] = result['gt']
        record.samples[sample].phased = result['phased']
        if 'PS' in record.header.formats:
            record.samples[sample]['PS'] = result['ps']
        statuses[sample] = result['status']
    return statuses
