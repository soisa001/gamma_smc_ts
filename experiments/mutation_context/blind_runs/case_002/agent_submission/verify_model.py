"""Independent string-based check of all context predictions and file checksums."""
from analyze import *
from itertools import product

summary=json.loads((ROOT/'summary.json').read_text())
rates={p['parameter']:p['rate'] for p in summary['parameters']}
table=np.load(ROOT/'final_context_probabilities.npy')
for i,tup in enumerate(product(BASES,repeat=7)):
    seq=''.join(tup);center=seq[3];p={base:0. for base in BASES}
    eq1='hot' if seq[0]==seq[5] else 'cold'
    eq2='hot' if seq[1]==seq[2] else 'cold'
    if center in 'AT':
        if center=='A':transition,transversion,other='G','C','T'
        else:transition,transversion,other='C','G','A'
        p[transition]=rates['at_transition_'+eq1]
        p[transversion]=rates['at_transversion_'+eq1]
        p[other]=rates['at_other_'+eq2]
    else:
        cp=(center=='C' and seq[4]=='G') or (center=='G' and seq[2]=='C')
        transition='T' if center=='C' else 'A'
        p[transition]=rates['cg_transition_'+('hot' if cp else 'cold')]
        for base in BASES:
            if base not in [center,transition]:p[base]=rates['cg_transversions']
    p[center]=1-sum(p.values())
    np.testing.assert_allclose(table[i],np.array([p[b] for b in BASES]),rtol=0,atol=1e-15)
for name,digest in summary['provenance']['input_sha256'].items():
    assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest
assert hashlib.sha256((ROOT/'predictions.tsv').read_bytes()).hexdigest()==summary['verification']['predictions_sha256']
print('PASS: all 16,384 context probability vectors match independently implemented string-based rules.')
print('PASS: all six original input hashes and the prediction hash match recorded provenance.')
