from mixture import *
def main():
    cc=np.load(ROOT/'counts11.npy',mmap_mode='r');k=11;mid=np.arange(4**k)//4**(k//2)%4
    r=192005/38400000;p=np.full((4**k,4),r/3);p[np.arange(4**k),mid]=1-r;out=[]
    for a in [10,30,100,300,1000,3000,10000]:
        q=(cc[0]+a*p)/(cc[0].sum(axis=1)[:,None]+a)
        row={'alpha':a,'tuning_ll':float(np.sum(cc[1]*np.log(q))),'assessment_ll':float(np.sum(cc[2]*np.log(q)))};out.append(row);print(row,flush=True)
    (ROOT/'validation11.json').write_text(json.dumps(out,indent=2))
if __name__=='__main__':main()
