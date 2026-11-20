from mixture import *
def main():
    p=fasta(ROOT/'calibration_parent.fasta');y=fasta(ROOT/'calibration_offspring.fasta');exp=json.loads((ROOT/'exploration.json').read_text());split={s['name']:s['split'] for s in exp['scaffolds']}
    shifts=list(range(-6,7));counts=np.zeros((13,3,4**9,4),dtype=np.uint32)
    for name,x in p.items():
        code=contexts(x,9)
        # Labels encode no change followed by cyclic alternate alleles, so windows need not include center.
        z=(y[name].astype(int)-x)%4
        for j,shift in enumerate(shifts):np.add.at(counts[j,split[name]],(np.roll(code,-shift),z),1)
    np.save(ROOT/'shift9_counts.npy',counts)
    out=[]
    for j,shift in enumerate(shifts):
        c=counts[j];w,ps,info=fit(c[0],1,steps=300);q=predict(c[0],w,ps)
        row={'shift':shift,'tuning_ll':float(np.sum(c[1]*np.log(q))),'assessment_ll':float(np.sum(c[2]*np.log(q))),**info};out.append(row);print(json.dumps(row),flush=True)
    (ROOT/'shift_validation.json').write_text(json.dumps(out,indent=2))
if __name__=='__main__':main()
