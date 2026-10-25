"""Fixed available-storm paired bootstrap, conditional on frozen predictions.

Seed 20261025 fixed before generating any revised interval. No fitting.
Input storm scores are four-lead case averages averaged within each storm.
"""
from pathlib import Path
import argparse, json, hashlib
import numpy as np
import pandas as pd

def main(data, out):
    out.mkdir(parents=True, exist_ok=True)
    st = pd.read_csv(data/'storm_scores.csv', float_precision='round_trip')
    old = pd.read_csv(data/'summary.csv', float_precision='round_trip')
    validation = [x for x in st.validation.unique() if x != 'LOSO']
    assert len(validation) == 1, validation
    st = st[st.validation == validation[0]]
    rows, sidmap, draws = [], {}, {}
    for center, n in [('ECMWF',45), ('NCEP',67)]:
        g = st[st.center == center]
        sids = sorted(g.sid.unique()); assert len(sids) == n
        sidmap[center] = sids
        draw = np.random.default_rng(20261025).integers(0,n,size=(5000,n))
        draws[center] = draw
        for (base,method,lead), z in g.groupby(['base','method','lead']):
            assert set(z.sid) == set(sids) and len(z) == n
            for metric in ['error_km','gain_vs_full_km','gain_vs_own_km','ATE_km','CTE_km']:
                if z[metric].isna().all() or (lead == '48_72_96_120_mean' and metric in ['ATE_km','CTE_km']): continue
                a = z.set_index('sid').loc[sids,metric].to_numpy()
                assert np.isfinite(a).all()
                mean = a.mean(); lo,hi = np.quantile(a[draw].mean(axis=1),[.025,.975])
                q = old[(old.center==center)&(old.validation==validation[0])&(old.base==base)&(old.method==method)&(old.lead==lead)&(old.metric==metric)]
                assert len(q)==1
                q=q.iloc[0]; assert abs(mean-q['mean']) < 1e-10
                rows.append(dict(center=center,validation=validation[0],base=base,method=method,lead=lead,metric=metric,mean=mean,ci_low=lo,ci_high=hi,old_ci_low=q.ci_low,old_ci_high=q.ci_high,cases=q.cases,storms=n))
    pd.DataFrame(rows).to_csv(out/'summary_fixed_cohort.csv', index=False,float_format='%.17g')
    np.savez_compressed(out/'paired_storm_draws.npz', **draws)
    (out/'storm_ids.json').write_text(json.dumps(sidmap,indent=2))
    report=dict(status='PASS',seed=20261025,replicates=5000,storms={k:len(v) for k,v in sidmap.items()},models_refit=0,point_means_unchanged=True,scope='Prior-year fixed available cohorts; conditional on saved predictions; paired across all methods and targets within center',input_sha256={n:hashlib.sha256((data/n).read_bytes()).hexdigest() for n in ['storm_scores.csv','summary.csv']})
    (out/'acceptance.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();main(a.data,a.out)
