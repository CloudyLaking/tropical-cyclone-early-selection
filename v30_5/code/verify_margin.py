"""Reconstruct paired circle-margin estimates."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--expected',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
sha=lambda x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
files=list(a.inputs.iterdir())+[a.expected];hashes={str(x):sha(x) for x in files if x.is_file()}
cases=pd.read_csv(a.inputs/'cases.csv.gz');expected=pd.read_csv(a.expected)
rows=cases[(cases.main806)&(cases.lead==120)]
assert len(rows)==4836 and not rows.duplicated(['code','protocol','group']).any()
checks=[]
for center,n in [('ECMWF',64),('NCEP',100)]:
 meta=json.loads((a.inputs/f'bootstrap_{center}.json').read_text());sids=meta['sids'];assert len(sids)==n
 draws=np.load(a.inputs/f'bootstrap_{center}.npy');assert draws.shape==(5000,n)
 for protocol in ['A','B']:
  z=rows[(rows.center==center)&(rows.protocol==protocol)]
  selected=z[z.group=='selected'].set_index('code').sort_index();random=z[z.group=='random'].set_index('code').sort_index()
  assert selected.index.equals(random.index)
  margin=(selected.r70-selected.distance).to_numpy()-(random.r70-random.distance).to_numpy()
  values=np.array([margin[selected.sid.to_numpy()==sid].mean() for sid in sids])
  mean=float(values.mean());lo,hi=np.percentile(values[draws].mean(axis=1),[2.5,97.5])
  ref=expected[(expected.center==center)&(expected.protocol==protocol)].iloc[0]
  err=float(np.max(np.abs([mean-ref['mean'],lo-ref.ci_low,hi-ref.ci_high])))
  assert err<1e-7,(center,protocol,err)
  checks.append(dict(center=center,protocol=protocol,mean=mean,ci_low=lo,ci_high=hi,max_error_km=err))
assert all(sha(f)==h for f,h in hashes.items())
a.out.mkdir(parents=True,exist_ok=False)
(a.out/'independent_margin_check.json').write_text(json.dumps(dict(status='PASS',input_sha256=hashes,source_sha256=sha(__file__),checks=checks,
 implementation='Direct per-case margin subtraction, explicit storm averaging, saved paired draws; no new random draws.'),indent=2))
print(json.dumps(dict(status='PASS',comparisons=len(checks),max_error_km=max(x['max_error_km'] for x in checks))))
