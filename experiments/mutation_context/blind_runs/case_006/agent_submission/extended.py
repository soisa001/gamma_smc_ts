from analyze import *
def main():
    p=fasta(ROOT/'calibration_parent.fasta');y=fasta(ROOT/'calibration_offspring.fasta');ex=json.loads((ROOT/'exploration.json').read_text());splits={s['name']:s['split'] for s in ex['scaffolds']}
    c=np.zeros((3,4**11,4),dtype=np.uint32)
    near=np.zeros(100,dtype=np.int64);adj=np.zeros(100,dtype=np.int64)
    for name,x in p.items():
        z=y[name];np.add.at(c[splits[name]],(contexts(x,11),z),1)
        m=x!=z
        for d in range(1,101):near[d-1]+=np.sum(m&np.roll(m,d))
        idx=np.flatnonzero(m);gaps=np.diff(np.r_[idx,idx[0]+len(x)])
        adj+=np.histogram(gaps,bins=np.r_[np.arange(1,100),100,np.inf])[0]
    np.save(ROOT/'counts11.npy',c)
    (ROOT/'spacing.json').write_text(json.dumps({'pairs_at_distances_1_to_100':near.tolist(),'gap_bins_1_to_99_and_100plus':adj.tolist()}))
    print('Finished 11-mer counts and event-spacing checks',near[:10])
if __name__=='__main__':main()
