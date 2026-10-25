"""Reconstruct frozen coordinate scores and paired statistical summaries."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import gzip,pickle,json,time,hashlib,traceback
import numpy as np
import pandas as pd

import argparse
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);args=p.parse_args()
R=Path(__file__).resolve().parents[1]
P=R/'data/protocol'
O=args.out
O.parent.mkdir(parents=True,exist_ok=True)
METRICS=['ATE','CTE','distance','r70','spread','in70','ellipse_in70','ellipse_radius','ellipse_rank0','ellipse_rank1']
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def score_batches(groups,truth,previous):
    x=np.asarray(groups,float)
    means=x.mean(axis=1)
    centered=x-means[:,None,:]
    member_radius=np.sqrt(np.sum(centered**2,axis=2))
    radius=np.quantile(member_radius,.7,axis=1,method='linear')
    delta=np.asarray(truth)-means
    distance=np.linalg.norm(delta,axis=1)
    direction=np.asarray(truth)-np.asarray(previous);direction/=np.linalg.norm(direction)
    right=np.array([direction[1],-direction[0]])
    # Sample-covariance eigenspace recovered through a centered-data SVD.
    _,singular,vh=np.linalg.svd(centered,full_matrices=False)
    eigen=singular**2/(x.shape[1]-1)
    tolerance=np.maximum(1e-12,1e-12*eigen.max(axis=1))
    active=eigen>tolerance[:,None]
    safe=np.where(active,eigen,1)
    components=np.einsum('gki,gji->gkj',centered,vh)
    md2=np.sum(components**2/safe[:,None,:]*active[:,None,:],axis=2)
    cutoff=np.quantile(md2,.7,axis=1,method='linear')
    target_components=np.einsum('gji,gi->gj',vh,delta)
    target_md2=np.sum(target_components**2/safe*active,axis=1)
    outside=np.sqrt(np.sum(target_components**2*(~active),axis=1))
    rank=active.sum(axis=1)
    ellipse=(target_md2<=cutoff)&(outside<=1e-8)
    ellipse[rank==0]=distance[rank==0]<=1e-8
    semiaxes=np.sqrt(np.maximum(cutoff,0)[:,None]*np.where(active,eigen,0))
    eqradius=np.sqrt(np.prod(semiaxes,axis=1))
    values=dict(ATE=-(delta@direction),CTE=-(delta@right),distance=distance,r70=radius,
                spread=np.sqrt(np.mean(member_radius**2,axis=1)),in70=(distance<=radius).astype(float),
                ellipse_in70=ellipse.astype(float),ellipse_radius=eqradius,
                ellipse_rank0=(rank==0).astype(float),ellipse_rank1=(rank==1).astype(float))
    margin=dict(circle=distance-radius,ellipse=target_md2-cutoff,outside=outside)
    return means,values,margin

def main():
    start=time.monotonic();O.mkdir(exist_ok=True)
    assert not (O/'acceptance.json').exists() and not (O/'failure.json').exists()
    with gzip.open(P/'frozen/cases.pkl.gz','rb') as f:records=pickle.load(f)
    with gzip.open(P/'frozen/random_lists.pkl.gz','rb') as f:random_lists=pickle.load(f)
    mainframe=pd.read_csv(P/'results/cases.csv.gz',float_precision='round_trip')
    mainsummary=pd.read_csv(P/'results/summary.csv',float_precision='round_trip')
    paths=[P/'frozen/cases.pkl.gz',P/'frozen/random_lists.pkl.gz',P/'results/cases.csv.gz',P/'results/summary.csv']
    hashes={str(p):sha(p) for p in paths}
    main_index=mainframe.set_index(['code','lead','protocol','group'])
    assert main_index.index.is_unique
    maincodes={p.stem for p in (R/'data/ranking/cases').glob('*.npz')}
    assert len(maincodes)==806
    rows=[];boundaries=[];degeneracy=[];errors=dict(mean_identity=0.,case_values=0.,summary=0.,interval=0.)
    keys=set()
    for count,rec in enumerate(records):
        code,lead=rec['code'],int(rec['lead']);keys.add((code,lead))
        assert bool(rec['main806'])==(code in maincodes)
        raw=np.asarray(rec['xy'],float);aligned=np.asarray(rec['translated_xy'],float);n=len(raw)
        selected=np.asarray(rec['selected'],int);random=np.asarray(random_lists[code],int)
        assert len(selected)==n//2 and len(np.unique(selected))==len(selected)
        assert random.shape==(128,len(selected)) and random.min()>=0 and random.max()<n
        assert np.all(np.diff(np.sort(random,axis=1),axis=1)>0)
        assert np.isfinite(raw).all() and np.isfinite(aligned).all()
        assert np.linalg.norm(np.asarray(rec['truth_xy'])-np.asarray(rec['truth_previous_xy']))>=1
        common=dict(code=code,center=rec['center'],sid=rec['sid'],lead=lead,main806=bool(rec['main806']))
        for group,indices in [('full',np.arange(n)[None,:]),('selected',selected[None,:]),('random',random)]:
            a=aligned[indices];b=raw[indices]+(a.mean(axis=1)-raw[indices].mean(axis=1))[:,None,:]
            meanerr=float(np.max(np.abs(a.mean(axis=1)-b.mean(axis=1))))
            errors['mean_identity']=max(errors['mean_identity'],meanerr)
            for protocol,positions in [('A',a),('B',b)]:
                means,values,margin=score_batches(positions,rec['truth_xy'],rec['truth_previous_xy'])
                reference=main_index.loc[code,lead,protocol,group]
                averaged={metric:float(values[metric].mean()) for metric in METRICS}
                for metric,value in averaged.items():
                    delta=abs(value-float(reference[metric]))
                    errors['case_values']=max(errors['case_values'],delta)
                    assert delta<=1e-7,(code,lead,protocol,group,metric,value,reference[metric])
                rows.append(dict(**common,protocol=protocol,group=group,**averaged))
                degeneracy.append(dict(**common,protocol=protocol,group=group,clouds=len(indices),rank0=int(values['ellipse_rank0'].sum()),rank1=int(values['ellipse_rank1'].sum())))
                for kind in ('circle','ellipse'):
                    close=np.flatnonzero(np.abs(margin[kind])<=1e-8)
                    for repeat in close:
                        boundaries.append(dict(**common,protocol=protocol,group=group,repeat=int(repeat),kind=kind,margin=float(margin[kind][repeat]),outside=float(margin['outside'][repeat]),covered=float(values['in70' if kind=='circle' else 'ellipse_in70'][repeat])))
        if count%500==0:print(json.dumps(dict(progress=count,total=len(records),elapsed_seconds=round(time.monotonic()-start,2),errors=errors)),flush=True)
    assert len(keys)==len(records) and len(records)==6885
    assert sum(bool(r['main806']) for r in records)==3224
    assert errors['mean_identity']<=1e-7
    frame=pd.DataFrame(rows);restored=[]
    for cohort,g in [('available',frame),('main806',frame[frame.main806])]:
        for (center,lead),cell in g.groupby(['center','lead']):
            fp=P/f'frozen/bootstrap_{cohort}_{center}_{lead}.npz';hashes[str(fp)]=sha(fp)
            with np.load(fp,allow_pickle=False) as bs:
                sids=bs['sids'].astype(str).tolist();draw=bs['draws']
            assert draw.shape==(5000,len(sids)) and sorted(cell.sid.unique())==sids
            if cohort=='main806':
                assert sids==json.loads((R/f'data/margin/bootstrap_{center}.json').read_text())['sids']
                assert np.array_equal(draw,np.load(R/f'data/margin/bootstrap_{center}.npy'))
            counts=np.stack([np.bincount(d,minlength=len(sids)) for d in draw])
            levels={}
            for (protocol,group),part in cell.groupby(['protocol','group']):
                levels[protocol+'_'+group]=np.stack([part.loc[part.sid==sid,METRICS].to_numpy().mean(axis=0) for sid in sids])
            comparisons=dict(levels)
            for protocol in ('A','B'):
                for baseline in ('full','random'):
                    comparisons[f'{protocol}_selected_minus_{baseline}']=levels[protocol+'_selected']-levels[protocol+'_'+baseline]
            for group in ('full','selected','random'):
                comparisons['A_minus_B_'+group]=levels['A_'+group]-levels['B_'+group]
            for baseline in ('full','random'):
                comparisons['A_minus_B_selected_minus_'+baseline]=comparisons['A_selected_minus_'+baseline]-comparisons['B_selected_minus_'+baseline]
            for comparison,stormvalues in comparisons.items():
                means=stormvalues.mean(axis=0);boot=counts@stormvalues/len(sids);ci=np.quantile(boot,[.025,.975],axis=0)
                for j,metric in enumerate(METRICS):
                    ref=mainsummary[(mainsummary.cohort==cohort)&(mainsummary.center==center)&(mainsummary.lead==lead)&(mainsummary.comparison==comparison)&(mainsummary.metric==metric)]
                    assert len(ref)==1,(cohort,center,lead,comparison,metric)
                    ref=ref.iloc[0]
                    errors['summary']=max(errors['summary'],abs(means[j]-float(ref['mean'])))
                    errors['interval']=max(errors['interval'],abs(ci[0,j]-float(ref.ci_low)),abs(ci[1,j]-float(ref.ci_high)))
                    restored.append(dict(cohort=cohort,center=center,lead=lead,comparison=comparison,metric=metric,mean=means[j],ci_low=ci[0,j],ci_high=ci[1,j]))
    assert errors['summary']<=1e-8 and errors['interval']<=1e-8,errors
    for path,digest in hashes.items():assert sha(path)==digest,path
    pd.DataFrame(restored).to_csv(O/'independent_summary.csv',index=False)
    pd.DataFrame(boundaries).to_csv(O/'boundary_cases.csv',index=False)
    deg=pd.DataFrame(degeneracy);descriptions=[]
    for cohort,g in [('available',deg),('main806',deg[deg.main806])]:
        q=g.groupby(['center','lead','protocol','group'])[['clouds','rank0','rank1']].sum().reset_index();q.insert(0,'cohort',cohort);descriptions.append(q)
    pd.concat(descriptions).to_csv(O/'degenerate_cloud_counts.csv',index=False)
    receipt=dict(status='PASS',case_leads=len(records),main_case_leads=3224,method_rows=len(frame),errors=errors,boundary_count=len(boundaries),
                 independent_ellipse='centered-data SVD',main_solver_imported=False,scope='Frozen derived coordinates, not raw archive reconstruction',source_hashes_unchanged=True,
                 script_sha256=sha(__file__),source_sha256=hashes,elapsed_seconds=time.monotonic()-start)
    (O/'acceptance.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in receipt.items() if k!='source_sha256'}),flush=True)

if __name__=='__main__':
    try:main()
    except BaseException:
        O.mkdir(exist_ok=True)
        (O/'failure.json').write_text(json.dumps(dict(status='FAIL',traceback=traceback.format_exc(),script_sha256=sha(__file__)),indent=2),encoding='utf-8')
        raise
