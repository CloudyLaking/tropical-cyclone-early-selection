"""Recheck frozen optimization certificates and reconstruct paired estimates.

For a convex set, z minimizes distance to zero iff z is feasible and
(v-z) dot z >= 0 for every vertex v. Monotone weights give the convex
hull of prefix (decreasing) or suffix (increasing) means. This verifies
every saved optimum without rerunning a search or drawing new permutations.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import argparse,json,time
import numpy as np
import pandas as pd

def main(out):
 root=Path(__file__).resolve().parents[1];data=root/'data';out.mkdir(parents=True,exist_ok=False)
 read=lambda p:pd.read_csv(data/p,float_precision='round_trip')
 rank=read('ranking/case_results.csv.gz').set_index(['code','lead']);paired=read('fig3/paired_cases.csv.gz').set_index(['code','lead'])
 ref_summary=read('ranking/summary.csv');contrast_summary=read('fig3/contrast_summary.csv');rows=[];pairs=[];certificates=0;started=time.monotonic()
 errors=dict(feasibility=0.,endpoint_km=0.,distance_km=0.,support_km2=0.,case_km=0.,summary_km=0.,ci_km=0.)
 def check(key,delta,tol):
  delta=float(delta);errors[key]=max(errors[key],delta);assert delta<=tol,(key,delta,tol)
 def segment(a,b):
  v=b-a;den=np.sum(v*v,axis=-1);u=np.divide(-np.sum(a*v,axis=-1),den,out=np.zeros_like(den),where=den>0)
  return np.linalg.norm(a+np.clip(u,0,1)[...,None]*v,axis=-1)
 for number,fp in enumerate(sorted((data/'ranking/cases').glob('*.npz'))):
  code=fp.stem
  with np.load(fp,allow_pickle=False) as z:
   n=len(z['ids']);assert np.array_equal(z['leads'],[48,72,96,120]);early=z['early_order'];random=z['random_orders']
   assert np.array_equal(np.sort(early),np.arange(n)) and np.array_equal(np.sort(random,axis=1),np.tile(np.arange(n),(128,1)))
   for j,lead in enumerate(z['leads']):
    p=z['error_xy'][j];late=z['late_orders'][j];assert np.array_equal(late,np.lexsort((z['ids'].astype(int),np.linalg.norm(p,axis=1))))
    check('endpoint_km',np.abs(z['translated_xy'][j]-z['target'][j]-p).max(),1e-7)
    vals={};h=z['H_weights'][j];hp=h@p
    check('feasibility',max(0.,-h.min(),abs(h.sum()-1)),1e-10)
    check('endpoint_km',np.linalg.norm(hp-z['H_points'][j]),1e-7)
    check('distance_km',abs(np.linalg.norm(hp)-z['H_distance'][j]),1e-7)
    check('support_km2',max(0.,-np.min((p-hp)@hp)),1e-5)
    vals['H']=float(np.linalg.norm(hp))
    for family,orders in [('early',early[None,:]),('late',late[None,:]),('random',random)]:
     for direction in [0,1]:
      w=z[family+'_weights'][j,:,direction];point=np.einsum('kn,nc->kc',w,p);ordered=np.take_along_axis(w,orders,axis=1)
      diff=np.diff(ordered,axis=1)*(1 if direction==0 else -1)
      check('feasibility',max(0.,-w.min(),np.abs(w.sum(axis=1)-1).max(),diff.max()),1e-10)
      check('endpoint_km',np.linalg.norm(point-z[family+'_points'][j,:,direction],axis=1).max(),1e-7)
      check('distance_km',np.abs(np.linalg.norm(point,axis=1)-z[family+'_distance'][j,:,direction]).max(),1e-7)
      op=p[orders[:,::-1] if direction else orders];vertices=op.cumsum(axis=1)/np.arange(1,n+1)[None,:,None]
      support=np.einsum('knc,kc->kn',vertices-point[:,None,:],point)
      check('support_km2',max(0.,-support.min()),1e-5);certificates+=len(orders)
     distance=z[family+'_distance'][j]
     assert np.array_equal(z[family+'_M_direction'][j],distance.argmin(axis=1))
     vals['D_'+family]=float(distance[:,0].mean());vals['I_'+family]=float(distance[:,1].mean());vals['M_'+family]=float(distance.min(axis=1).mean())
     vals['I_better_'+family]=float(np.mean(distance[:,1]<distance[:,0]-1e-9))
    for direction in ['D','M']:
     for label,left,right in [('early_minus_late','early','late'),('late_minus_H','late',None),('early_minus_random','early','random'),('early_minus_H','early',None),('random_minus_H','random',None)]:
      v=vals[direction+'_'+left]-(vals[direction+'_'+right] if right else vals['H']);name=direction+'_'+label
      vals[name]=v;vals[name+'_positive_1km']=float(v>1);vals[name+'_negative_1km']=float(v < -1)
    ref=rank.loc[code,lead]
    # I_better indicators are preserved with the source tolerance, not redefined here.
    for k in list(vals):
     if k.startswith('I_better_'):del vals[k];continue
     check('case_km',abs(vals[k]-ref[k]),1e-7)
    meta=dict(code=code,lead=int(lead),center=ref.center,sid=ref.sid);rows.append(dict(**meta,**vals))
    op=p[early];sums=op.cumsum(axis=0);k=np.arange(1,n);a=sums[:-1]/k[:,None];b=(sums[-1]-sums[:-1])/(n-k)[:,None];ds=segment(a,b)
    actual=float(np.linalg.norm(op[:n//2].mean(axis=0)));share=float(ds[n//2-1]);split=float(ds.min())
    pv=dict(actual=actual,fixed_half_share=share,two_level=split,half_minus_share=actual-share,share_minus_split=share-split,split_minus_levels=split-vals['M_early'],early_minus_verifying=vals['M_early']-vals['M_late'],verifying_minus_free=vals['M_late']-vals['H'])
    for k,v in pv.items():check('case_km',abs(v-paired.loc[code,lead][k]),1e-7)
    pairs.append(dict(**meta,**pv))
  if number%100==0:print('Checked',number+1,'/806',flush=True)
 assert len(rows)==3224 and certificates==838240
 frame=pd.DataFrame(rows);pairframe=pd.DataFrame(pairs);restored=[];restored_pairs=[]
 for center,g in frame.groupby('center'):
  sids=json.loads((data/f'margin/bootstrap_{center}.json').read_text())['sids'];draws=np.load(data/f'margin/bootstrap_{center}.npy');assert draws.shape==(5000,len(sids))
  metrics=[c for c in frame if c not in ['code','lead','center','sid']]
  cells=list(g.groupby('lead'))+[('48-120h mean',g.groupby(['code','sid'],as_index=False)[metrics].mean())]
  for lead,part in cells:
   storm=part.groupby('sid')[metrics].mean().loc[sids].to_numpy();mean=storm.mean(axis=0);ci=np.quantile(storm[draws].mean(axis=1),[.025,.975],axis=0)
   for k,metric in enumerate(metrics):
    ref=ref_summary[(ref_summary.center==center)&(ref_summary.lead.astype(str)==str(lead))&(ref_summary.metric==metric)].iloc[0]
    check('summary_km',abs(mean[k]-ref['mean']),1e-8);check('ci_km',max(abs(ci[0,k]-ref.ci_low),abs(ci[1,k]-ref.ci_high)),1e-8)
    restored.append(dict(center=center,lead=lead,metric=metric,mean=mean[k],ci_low=ci[0,k],ci_high=ci[1,k]))
  for metric in ['half_minus_share','share_minus_split','split_minus_levels','early_minus_verifying','verifying_minus_free']:
   v=pairframe[(pairframe.center==center)&(pairframe.lead==120)].groupby('sid')[metric].mean().loc[sids].to_numpy();mean=v.mean();ci=np.quantile(v[draws].mean(axis=1),[.025,.975]);ref=contrast_summary[(contrast_summary.center==center)&(contrast_summary.contrast==metric)].iloc[0]
   check('summary_km',abs(mean-ref['mean']),1e-8);check('ci_km',max(abs(ci[0]-ref.ci_low),abs(ci[1]-ref.ci_high)),1e-8)
   restored_pairs.append(dict(center=center,contrast=metric,mean=mean,ci_low=ci[0],ci_high=ci[1]))
 pd.DataFrame(restored).to_csv(out/'ranking_summary.csv',index=False);pd.DataFrame(restored_pairs).to_csv(out/'paired_summary.csv',index=False)
 receipt=dict(status='PASS',cases=806,case_leads=3224,ordered_certificates=certificates,errors=errors,seconds=time.monotonic()-started,scope='Frozen derived inputs. All saved monotone and free optima checked by convex projection inequalities; all two-level segments recomputed. No raw-source reconstruction, new permutations or bootstrap draws.')
 (out/'acceptance.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);main(p.parse_args().out)
