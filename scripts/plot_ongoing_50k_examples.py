#!/usr/bin/env python3
"""Prepare and plot 50-kya-to-present examples at the selected site.

Fixed scientific seeds; first accepted s=.01 proposal. The optional s=.05
completed sweep replays proposal 0 to present, only when fixation is allowed.
Phases: simulate-fixed, prepare, report. Report uses committed caches only.
"""
import os
os.environ['MPLBACKEND'] = 'Agg'
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import re
import shutil
import subprocess
import time
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
import msprime
import numpy as np
import pandas as pd
from scipy.special import gammainc
from scipy.stats import gamma
import tskit
import zstandard
from gamma_smc_aou.fresh_power import ordered_nodes
from plot_fitted_gamma_examples import read_exact
from plot_selection_episode_series import sha, save_json
from report_selection_arrays import local_path
import run_selection_arrays as arrays

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'docs/results/selection_arrays_eas_20261105_4h'
NEUTRAL = ROOT / 'docs/results/selected_site_posteriors_20261106'
OLD_RAW = ROOT.parent / 'sim/selection_arrays_20261104'
NEW_RAW = ROOT.parent / 'sim/ongoing_50k_fixed_example_20261109'
FIXED_UNIT = '50k/introgression/s0.05/selected_00000'


def simulate_fixed():
    """Replay the existing first proposal, retaining its exact trajectory to fixation."""
    source = OLD_RAW / FIXED_UNIT
    original = json.loads((source/'rejected.json').read_text())
    assert original['rejection']=='fixation' and original['replicate']==0
    cfg = json.loads((OLD_RAW/'manifest.json').read_text())['contract']['configs']['50k']
    for key,path in [('slim_sha256',cfg['slim']),('decoder_sha256',cfg['decoder']),
                     ('score_source_sha256',arrays.base.__file__)]:
        assert sha(Path(path)) == cfg[key]
    destination = NEW_RAW / FIXED_UNIT
    destination.mkdir(parents=True,exist_ok=True)
    contract = dict(config=cfg, original_receipt_sha256=sha(source/'rejected.json'),
        original_model_sha256=sha(source/'model.slim'), original_onset_sha256=original['onset_sha256'],
        conditioning='Allow fixation; require focal allele present at end. Replay first proposal, no ascertainment by inferred signal.',
        change='Remove fixation early-stop from original SLiM model; retain loss stop, seed, initial state and every other model statement.')
    manifest = NEW_RAW/'manifest.json'
    if manifest.exists():assert json.loads(manifest.read_text())==contract
    else:save_json(manifest,contract)
    for name in ('pairs.tsv','positions.txt'):
        if (NEW_RAW/name).exists():assert sha(NEW_RAW/name)==sha(OLD_RAW/name)
        else:shutil.copy2(OLD_RAW/name,NEW_RAW/name)
    if (destination/'complete.json').exists():
        record=json.loads((destination/'complete.json').read_text())
        for name,key in [('trajectory.tsv','trajectory_sha256'),('observed.trees','observed_sha256'),
                         ('posterior.zst','posterior_sha256'),('scores.npz','scores_sha256')]:
            assert sha(destination/name)==record[key]
        return
    label=cfg['seed_namespace']+original['unit']
    receipt=destination/'simulation.json'
    if receipt.exists():
        record=json.loads(receipt.read_text())
        assert sha(destination/'observed.trees')==record['observed_sha256']
    else:
        assert sha(source/'onset.trees')==original['onset_sha256']
        shutil.copy2(source/'onset.trees',destination/'onset.trees')
        model=(source/'model.slim').read_text()
        old_stop='if (af==0.0 | af==1.0) sim.simulationFinished();'
        assert model.count(old_stop)==1
        model=model.replace(old_stop,'if (af==0.0) sim.simulationFinished();')
        model=model.replace(str(source),str(destination))
        assert str(source) not in model
        (destination/'model.slim').write_text(model)
        start=time.monotonic()
        result=subprocess.run([cfg['slim'],'-s',str(arrays.base.seed(label+':slim')),str(destination/'model.slim')],capture_output=True,text=True)
        (destination/'slim.stdout.log').write_text(result.stdout)
        (destination/'slim.stderr.log').write_text(result.stderr)
        result.check_returncode()
        raw=np.array([list(map(float,line.split()[1:])) for line in result.stdout.splitlines() if line.startswith('TRAJECTORY ')])
        np.testing.assert_array_equal(raw[:,0],np.arange(1,2002))
        data=np.c_[(2001-raw[:,0])*25,raw[:,4]/(2*raw[:,2]),raw[:,2:]]
        np.testing.assert_allclose(data[:,3],.05)
        prior=np.loadtxt(source/'trajectory.tsv',skiprows=1)
        np.testing.assert_allclose(data[:len(prior)],prior,rtol=0,atol=1e-14)
        assert data[-1,1]==1
        np.savetxt(destination/'trajectory.tsv',data,delimiter='\t',header='years_ago\taf\tdiploids\tselection_coefficient\tfocal_copies',comments='')
        ts=arrays.base._sample_diploids(tskit.load(destination/'forward.trees'),cfg['diploids'],arrays.base.seed(label+':sample'))
        tables=ts.dump_tables();tables.sites.clear();tables.mutations.clear();tables.time_units='generations'
        observed=msprime.sim_mutations(tables.tree_sequence(),rate=cfg['mutation_rate'],random_seed=arrays.base.seed(label+':mutations'),model=msprime.JC69())
        observed.dump(destination/'observed.trees')
        record=dict(original,accepted=True,rejection=None,stop_years_ago=0,trajectory_rows=len(data),
            conditioning=contract['conditioning'],source_directory=str(destination),
            original_receipt_sha256=contract['original_receipt_sha256'],config_sha256=sha(manifest),
            simulation_seconds=time.monotonic()-start,trajectory_sha256=sha(destination/'trajectory.tsv'),
            observed_sha256=sha(destination/'observed.trees'),first_fixation_kya=float(data[data[:,1]==1,0][0]/1000))
        save_json(receipt,record)
    native_lib=Path(cfg['slim']).parent.parent/'lib'
    os.environ['LD_LIBRARY_PATH']=str(native_lib)+(':'+os.environ['LD_LIBRARY_PATH'] if os.environ.get('LD_LIBRARY_PATH') else '')
    linkage=subprocess.run(['ldd',cfg['decoder']],capture_output=True,text=True,check=True)
    assert 'not found' not in linkage.stdout
    save_json(NEW_RAW/'execution_environment.json',dict(native_lib=str(native_lib),decoder_linkage=linkage.stdout,
        script_sha256=sha(Path(__file__)),phase_seeds={phase:arrays.base.seed(label+':'+phase) for phase in ('slim','sample','mutations')}))
    record.update(arrays.base.score(cfg,destination,NEW_RAW))
    save_json(destination/'complete.json',record)
    print(f"Completed s=.05 replay; fixed at {record['first_fixation_kya']} kya.",flush=True)


def extract(row):
    directory=local_path(row['source_directory'])
    for name,key in [('posterior.zst','posterior_sha256'),('scores.npz','scores_sha256'),
                     ('observed.trees','observed_sha256'),('trajectory.tsv','trajectory_sha256')]:
        assert sha(directory/name)==row[key]
    meta=json.loads((directory/'posterior.zst.meta').read_text())
    np.testing.assert_array_equal(meta['output_positions'],np.arange(100)*10000)
    assert meta['num_pairs']==1000 and meta['sequence_length']==100
    pairs=np.array(meta['pairs']);chunk=meta['chunk_size'];size=100*chunk*4;parts=[]
    with (directory/'posterior.zst').open('rb') as handle,zstandard.ZstdDecompressor().stream_reader(handle) as reader:
        for first in range(0,1000,chunk):
            n=min(chunk,1000-first)
            a=np.frombuffer(read_exact(reader,size),np.float32).reshape(100,chunk)[50,:n].astype(float)
            b=np.frombuffer(read_exact(reader,size),np.float32).reshape(100,chunk)[50,:n].astype(float)
            assert np.isfinite(a).all() and np.isfinite(b).all() and a.min()>0 and b.min()>0
            parts.append((a,b))
        assert not reader.read(1)
    a=np.concatenate([p[0] for p in parts])
    scale=meta['scaled_mutation_rate']/(2*1.29e-8)*.025/np.concatenate([p[1] for p in parts])
    ts=tskit.load(directory/'observed.trees');nodes=ordered_nodes(ts);tree=ts.at(500000)
    for i in range(200):assert meta['sample_names'][str(i)]==f'tsk_{i//2}.{i%2}'
    truth=np.array([tree.tmrca(int(nodes[x]),int(nodes[y]))*.025 for x,y in pairs])
    probabilities=gammainc(a[:,None],np.arange(1,51)/scale[:,None])
    with np.load(directory/'scores.npz') as data:cdfs=data['cdfs'][:,50].copy()
    for reconstructed,saved in zip([(truth[:,None]<np.arange(1,51)).mean(0),(probabilities>.5).mean(0),probabilities.mean(0)],cdfs):
        np.testing.assert_allclose(reconstructed,saved,rtol=1e-12,atol=1e-14)
    trajectory=np.loadtxt(directory/'trajectory.tsv',skiprows=1)
    np.testing.assert_array_equal(trajectory[:,0],np.arange(50000,-1,-25))
    np.testing.assert_allclose(trajectory[:,3],float(row['selection']))
    return a,scale,truth,pairs,cdfs,trajectory


def prepare(out,workers,include_fixed):
    inventory=json.loads((BASE/'output_hashes.json').read_text())
    assert sha(BASE/'completed_units.tsv')==inventory['completed_units.tsv']
    rows=pd.read_csv(BASE/'completed_units.tsv',sep='\t')
    chosen=rows[(rows.regime=='50k')&(rows.origin=='introgression')&(rows.cohort=='selected')&(rows.selection==.01)].sort_values('replicate').iloc[0].to_dict()
    assert chosen['replicate']==1 and 0<chosen['final_frequency']<1
    selected=[chosen]
    if include_fixed:
        record=json.loads((NEW_RAW/FIXED_UNIT/'complete.json').read_text())
        selected.append(record)
    neutral_meta=json.loads((NEUTRAL/'manifest.json').read_text())
    neutral_af=json.loads((NEUTRAL/'neutral_af_manifest.json').read_text())
    assert sha(NEUTRAL/'parameters.npz')==neutral_meta['cache_sha256']
    assert sha(NEUTRAL/'neutral_af.npz')==neutral_af['cache_sha256']
    contract=dict(selected_ids=[r['id'] for r in selected],selected_scores=[r['scores_sha256'] for r in selected],
        neutral_parameters_sha256=neutral_meta['cache_sha256'],neutral_af_sha256=neutral_af['cache_sha256'],
        position_bp=500000,selection_window_kya=[50,0],selection=[r['selection'] for r in selected],
        choice='s=.01: first surviving segregating proposal; s=.05 if included: first proposal replayed allowing fixation, unchosen by inferred signal.',
        schema='ongoing-50k-examples/v1')
    if (out/'manifest.json').exists():
        saved=json.loads((out/'manifest.json').read_text())
        assert saved['contract']==contract and sha(out/'parameters.npz')==saved['cache_sha256']
        return
    with ThreadPoolExecutor(max_workers=workers) as pool:parts=list(pool.map(extract,selected))
    with np.load(NEUTRAL/'parameters.npz') as old:
        for p in parts:np.testing.assert_array_equal(p[3],old['pairs'])
        cache={name:np.concatenate([np.stack([p[index] for p in parts]),old[name][3:]]) for name,index in [('alpha',0),('scale_kyr',1),('true_kya',2),('cdfs',4)]}
        ids=np.concatenate([np.array(contract['selected_ids']),old['ids'][3:]])
    with np.load(NEUTRAL/'neutral_af.npz') as old:
        np.testing.assert_array_equal(ids[len(selected):],old['ids'])
        neutral_trajectories=old['trajectories'].copy()
    np.savez_compressed(out/'parameters.npz',**cache,ids=ids,trajectories=np.stack([p[5] for p in parts]),neutral_trajectories=neutral_trajectories,pairs=parts[0][3])
    # Convert pandas missing values to JSON null, retaining complete receipt provenance.
    selected=json.loads(pd.DataFrame(selected).to_json(orient='records',double_precision=15))
    save_json(out/'manifest.json',dict(contract=contract,selected=selected,neutral_ids=ids[len(selected):].tolist(),
        cache_sha256=sha(out/'parameters.npz'),configuration=json.loads((BASE/'manifest.json').read_text())['contract']['configs']['50k'],
        validation='Raw hashes, exact 500-kb position, all 50 truth/median/mass cutoffs, identical 1000 pairs and all-generation selection coefficients checked.'))
    print(f'Prepared {len(selected)} selected examples and 100 neutrals.',flush=True)


def report(out,workers):
    meta=json.loads((out/'manifest.json').read_text());n=len(meta['selected'])
    assert sha(out/'parameters.npz')==meta['cache_sha256']
    with np.load(out/'parameters.npz') as data:
        a,scales,cdfs,trajectories,neutral_af=[data[k].copy() for k in ('alpha','scale_kyr','cdfs','trajectories','neutral_trajectories')]
    matplotlib.rcParams.update({'font.size':15,'axes.labelsize':16,'axes.titlesize':16,'xtick.labelsize':13,'ytick.labelsize':13,'pdf.fonttype':42,'ps.fonttype':42})
    times=np.linspace(.05,200,1000)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        densities=np.stack(list(pool.map(lambda i:gamma.pdf(times[:,None],a=a[i],scale=scales[i]).mean(1),range(n+100))))
    cumulative=np.concatenate([np.zeros((n+100,3,1)),cdfs],axis=2);bins=np.diff(cumulative,axis=2);edges=np.arange(51)
    colors=['#D55E00','#CC79A7']
    def background(ax,x,neutral,selected,color):
        for y in neutral:ax.plot(x,y,color='.6',lw=.7,alpha=.4)
        ax.plot(x,neutral.mean(0),color='#0072B2',lw=2.5)
        ax.plot(x,selected,color=color,lw=3)
        ax.yaxis.set_major_formatter(PercentFormatter(1));ax.grid(axis='y',alpha=.2);ax.set_ylim(bottom=0)
    def save(fig,stem):
        for extension in ('png','pdf'):fig.savefig(out/f'{stem}.{extension}',dpi=180)
        plt.close(fig)
    summary=[]
    for i,record in enumerate(meta['selected']):
        s=float(record['selection']);stem=f's{s:g}'
        status='segregating at present' if trajectories[i,-1,1]<1 else 'fixed at present'
        title=f'50 kya - present · s={s:g} · selected site (500 kb)'
        handles=[Line2D([],[],color='.6',lw=1,label='100 neutrals'),Line2D([],[],color='#0072B2',lw=2.5,label='Neutral mean'),Line2D([],[],color=colors[i],lw=3,label=f's={s:g}')]
        fig,axes=plt.subplots(2,1,figsize=(11,8.5))
        background(axes[0],times,densities[n:],densities[i],colors[i]);axes[0].set(xlim=(0,200),xlabel='Pairwise TMRCA (kya)',ylabel='Posterior density\n(% per kyr)')
        background(axes[1],edges,cumulative[n:,2],cumulative[i,2],colors[i]);axes[1].set(xlim=(0,50),xlabel='TMRCA cutoff (kya)',ylabel='Cumulative posterior\ncoalescence fraction')
        fig.suptitle(title,y=.98,fontsize=18);fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.94),ncol=3,frameon=False,fontsize=12)
        fig.text(.5,.025,f'EAS + 1% introgression · 1,000 pairs/simulation · {status} · no renormalization',ha='center',fontsize=11)
        fig.subplots_adjust(top=.84,bottom=.12,left=.16,right=.97,hspace=.6);save(fig,stem+'_posterior')
        fig,axes=plt.subplots(3,1,figsize=(11,8.5),sharex=True)
        for j,ax in enumerate(axes):
            background(ax,edges[:-1]+.5,bins[n:,j],bins[i,j],colors[i]);ax.set_ylabel('Pair fraction');ax.set_xlim(0,50)
            ax.set_title(['True TMRCA','Posterior median','Posterior-weighted mass'][j],loc='left',fontsize=14)
        axes[-1].set_xlabel('Pairwise TMRCA: 1-kyr intervals (kya)')
        fig.suptitle(title,y=.98,fontsize=18);fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.94),ncol=3,frameon=False,fontsize=12)
        fig.text(.5,.025,f'Same 1,000 pairs · 50 held-out + 50 calibration neutrals · {status}',ha='center',fontsize=11)
        fig.subplots_adjust(top=.82,bottom=.12,left=.13,right=.97,hspace=.57);save(fig,stem+'_intervals')
        summary.append(dict(selection=s,id=record['id'],present_af=float(trajectories[i,-1,1]),status=status,
            below_50_kya={method:dict(selected=float(cdfs[i,j,-1]),neutral_mean=float(cdfs[n:,j,-1].mean()),neutral_ge_selected=int(np.sum(cdfs[n:,j,-1]>=cdfs[i,j,-1]))) for j,method in enumerate(('truth','median','mass'))},
            posterior_tail_above_200=float(gamma.sf(200,a=a[i],scale=scales[i]).mean())))
    fig,ax=plt.subplots(figsize=(11,8.5));time_af=neutral_af[0,:,0]/1000
    for trajectory in neutral_af:ax.plot(time_af,trajectory[:,1],color='.6',lw=.8,alpha=.45)
    handles=[Line2D([],[],color='.6',lw=1,label='100 neutrals')]
    handles+=ax.plot(time_af,neutral_af[:,:,1].mean(0),color='#0072B2',lw=2.5,label='Neutral mean')
    for i,record in enumerate(meta['selected']):handles+=ax.plot(trajectories[i,:,0]/1000,trajectories[i,:,1],color=colors[i],lw=3,label=f"s={float(record['selection']):g}")
    ax.set(xlim=(50,0),ylim=(0,1.02),xlabel='Time before present (kya)',ylabel='Allele frequency');ax.yaxis.set_major_formatter(PercentFormatter(1));ax.grid(axis='y',alpha=.2)
    fig.suptitle('50 kya - present',y=.96,fontsize=21);fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.91),ncol=2,frameon=False,fontsize=16)
    fig.text(.5,.055,'EAS + 1% introgression at 50 kya · focal site: 500 kb',ha='center',fontsize=13)
    fig.text(.5,.025,'s=0.01: segregation required; s=0.05: fixation allowed. Neutrals unconditioned.' if n==2 else 's=0.01: present-day segregation required. Neutrals unconditioned.',ha='center',fontsize=11)
    fig.subplots_adjust(top=.78,bottom=.15,left=.13,right=.97);save(fig,'allele_frequencies')
    for path in out.glob('*.pdf'):
        contents=path.read_bytes();assert b'/FontFile2' in contents and not re.search(rb'/Subtype\s*/Image',contents)
    save_json(out/'summary.json',dict(examples=summary,script_sha256=sha(Path(__file__)),caveat='Illustrative examples, not power estimates. Fixed-site comparisons use 50 held-out and 50 calibration neutrals.'))
    (out/'.gitattributes').write_text('* -text\n')
    save_json(out/'output_hashes.json',{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='output_hashes.json'})
    print(json.dumps(summary,indent=2),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase',choices=['simulate-fixed','prepare','report','all'])
    parser.add_argument('--include-fixed',action='store_true')
    parser.add_argument('--out',type=Path,default=ROOT/'docs/results/ongoing_50k_examples_20261109')
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args();assert 1<=args.workers<=24
    args.out.mkdir(parents=True,exist_ok=True)
    if args.phase=='simulate-fixed':
        if not args.include_fixed:parser.error('Explicit --include-fixed required for completed sweep.')
        simulate_fixed()
    if args.phase in ('prepare','all'):prepare(args.out,args.workers,args.include_fixed)
    if args.phase in ('report','all'):report(args.out,args.workers)


if __name__=='__main__':main()
