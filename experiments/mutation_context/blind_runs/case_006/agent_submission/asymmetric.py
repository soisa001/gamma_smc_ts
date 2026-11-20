from mixture import *
def main():
    c11=np.load(ROOT/'counts11.npy',mmap_mode='r');out=[]
    for left,right in [(3,4),(4,3),(4,4),(4,5),(5,4),(5,5)]:
        k=left+right+1
        cc=c11.reshape(3,4**(5-left),4**k,4**(5-right),4).sum(axis=(1,3))
        mid=np.arange(4**k)//4**right%4
        c=np.take_along_axis(cc,(mid[:,None]+np.arange(4))[None,:,:]%4,axis=2)
        for h in [1,3]:
            w,ps,info=fit(c[0],h,steps=300);p=predict(c[0],w,ps)
            row={'left':left,'right':right,'h':h,'tuning_ll':float(np.sum(c[1]*np.log(p))),'assessment_ll':float(np.sum(c[2]*np.log(p))),**info};out.append(row);print(json.dumps(row),flush=True)
    (ROOT/'asymmetric_validation.json').write_text(json.dumps(out,indent=2))
if __name__=='__main__':main()
