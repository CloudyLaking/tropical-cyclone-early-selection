"""Verify frozen inputs, reconstruct principal estimates, and render main figures."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys,time
ROOT=Path(__file__).resolve().parent
def main():
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true',help='Check file hashes only; standard library only')
 p.add_argument('--checks',action='store_true',help='Reconstruct protocol, weight and margin estimates before plotting')
 p.add_argument('--out',type=Path,default=ROOT/'results');p.add_argument('--font-dir',type=Path,help='Optional directory with licensed arial.ttf and arialbd.ttf; otherwise DejaVu Sans')
 a=p.parse_args();manifest=json.loads((ROOT/'MANIFEST.json').read_text());bad=[f for f,h in manifest.items() if not (ROOT/f).is_file() or hashlib.sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 assert not bad,bad;print('Input/source SHA-256 verified:',len(manifest),'files',flush=True)
 if a.verify:return
 out=a.out.resolve();assert not out.exists(),'Choose a new --out directory; existing outputs are preserved';out.mkdir(parents=True)
 env={**os.environ,'REPRO_OUTPUT':str(out),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
 if a.font_dir:env['REPRO_FONT_DIR']=str(a.font_dir.resolve())
 commands=[]
 if a.checks:
  commands.extend([['code/verify_weights.py','--out',str(out/'weights')],['code/verify_protocol.py','--out',str(out/'protocol')],['code/verify_margin.py','--inputs','data/margin','--expected','data/fig4/margin_summary.csv','--out',str(out/'margin')]])
 commands.extend([['code/plot_schematic.py'],['code/plot_later_effects.py','--data','data/fig2','--out',str(out/'fig2')],['code/plot_constraints.py','--data','data/fig3','--out',str(out/'fig3')],['code/plot_ranges_landscape.py','--data','data/fig4','--out',str(out/'fig4')]])
 start=time.monotonic()
 for cmd in commands:
  print('Running',cmd[0],flush=True);subprocess.run([sys.executable]+cmd,cwd=ROOT,env=env,check=True)
 (out/'run_receipt.json').write_text(json.dumps(dict(status='PASS',checks=a.checks,commands=commands,source_manifest_sha256=hashlib.sha256((ROOT/'MANIFEST.json').read_bytes()).hexdigest(),seconds=time.monotonic()-start),indent=2))
 print('Complete:',out,flush=True)
if __name__=='__main__':main()
