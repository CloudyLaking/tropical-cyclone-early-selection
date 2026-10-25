"""Derived-data reconstruction: horizontal display of frozen clouds and paired diagnostics."""
from pathlib import Path
import argparse,json,hashlib
from plot_fonts import configure_font
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from audit_panel_alignment import require_matplotlib_panel_alignment
SHA=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
T,O,P='#0095A6','#E98B2A','#8B5FBF'
def main(data,out,fonts):
 out.mkdir(parents=True,exist_ok=False)
 hashes={p.name:SHA(p) for p in data.iterdir() if p.is_file()}
 FONT=configure_font(fonts)
 plt.rcParams.update({'font.family':FONT,'font.size':8.5,'axes.labelsize':8.5,'axes.titlesize':8.5,'xtick.labelsize':8.5,'ytick.labelsize':8.5,'legend.fontsize':8.5,'axes.linewidth':.65,'pdf.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white','axes.facecolor':'white','axes.spines.top':True,'axes.spines.right':True,'legend.frameon':False})
 read=lambda n:pd.read_csv(data/n)
 m=read('source_cloud_members.csv');b=read('source_region_boundaries.csv');g=read('source_region_geometry.csv');dr=read('source_center_radius_changes.csv');cov=read('source_paired_coverage.csv')
 assert len(m)==70 and len(cov)==16 and len(dr)==8
 assert json.loads((data/'margin_acceptance.json').read_text())['status']=='PASS'
 W,H=183,120
 fig=plt.figure(figsize=(W/25.4,H/25.4))
 axes=[fig.add_axes([x/W,y/H,w/W,31/H]) for y in (75,27) for x,w in [(14,44),(80,37),(143,35)]]
 def title(ax,l,t):
  ax.text(0,1.07,l,transform=ax.transAxes,fontsize=8.5,fontweight='bold',va='bottom')
  ax.text(.10,1.07,t,transform=ax.transAxes,fontsize=8.5,va='bottom')
 lim=600
 for i,center in enumerate(('ECMWF','NCEP')):
  cloud,change,coverage=axes[i*3:i*3+3]
  for protocol,color,marker in [('MW',T,'o'),('CT',O,'^')]:
   z=m[(m.center==center)&(m.protocol==protocol)]
   cloud.scatter(z.relative_east_km,z.relative_north_km,s=8,marker=marker,color=color,alpha=.55,linewidths=0,zorder=3)
   for shape,ls in [('circle','-'),('ellipse','--')]:
    q=b[(b.center==center)&(b.protocol==protocol)&(b['shape']==shape)].sort_values('node')
    cloud.plot(q.east_km,q.north_km,color=color,ls=ls,lw=.8,zorder=2)
  z=m[(m.center==center)&(m.protocol=='MW')]
  truth=z[['truth_relative_east_km','truth_relative_north_km']].iloc[0].to_numpy()
  radius=float(g[(g.center==center)&(g.protocol=='MW')].iloc[0].r70_km)
  cloud.plot([0,truth[0]],[0,truth[1]],color='#555555',lw=.7,zorder=4)
  cloud.annotate('d',xy=(truth[0]*.58,truth[1]*.58),xytext=(400,400),fontsize=8.5,color='#333333',arrowprops=dict(arrowstyle='-',color='#777777',lw=.5))
  cloud.plot([0,radius],[0,0],color='#555555',lw=.7,zorder=4)
  cloud.annotate('r',xy=(radius*.7,0),xytext=(400,-350),fontsize=8.5,color='#333333',arrowprops=dict(arrowstyle='-',color='#777777',lw=.5))
  cloud.plot(0,0,marker='s',ms=4,mfc='white',mec='black',mew=.8,zorder=6)
  cloud.plot(*truth,marker='+',ms=6,color='black',mew=1,zorder=6)
  cloud.set(xlim=(-lim*44/31,lim*44/31),ylim=(-lim,lim),xticks=[-400,0,400],yticks=[-400,0,400],xlabel='East of selected mean (km)',ylabel='North of selected mean (km)')
  cloud.set_aspect('equal',adjustable='box')
  cloud.text(.03,.97,'Circle: covered if d ≤ r',transform=cloud.transAxes,ha='left',va='top',fontsize=8.5)
  title(cloud,'ad'[i],center+' · same mean')
  change.axvline(0,ymax=.80,color='#777777',ls='--',lw=.65)
  for j,protocol in enumerate(('A','B')):
   for metric,color,offset,marker in [('distance',T,.12,'s'),('r70',O,-.12,'o')]:
    q=dr[(dr.center==center)&(dr.protocol==protocol)&(dr.metric==metric)].iloc[0]
    change.errorbar(q['mean'],1-j+offset,xerr=[[q['mean']-q.ci_low],[q.ci_high-q['mean']]],fmt=marker,ms=3.6,color=color,capsize=1.8,lw=.8)
  change.set(yticks=[1,0],yticklabels=['MW','CT'],ylim=(-.5,1.5),xlim=(-60,20),xticks=[-60,-30,0],xlabel='Selected − random (km)')
  title(change,'be'[i],'Center / radius')
  coverage.axvline(0,ymax=.82 if center=='NCEP' else 1,color='#777777',ls='--',lw=.65)
  configs=[('A','in70','MW circle'),('B','in70','CT circle'),('A','ellipse_in70','MW ellipse'),('B','ellipse_in70','CT ellipse')]
  cohorts=[('main806',T,.13,'o'),('available',P,-.13,'D')] if center=='ECMWF' else [('main806',T,0,'o')]
  for cohort,color,offset,marker in cohorts:
   for j,(protocol,metric,label) in enumerate(configs):
    q=cov[(cov.center==center)&(cov.cohort==cohort)&(cov.comparison==protocol+'_selected_minus_random')&(cov.metric==metric)].iloc[0]
    mean,lo,hi=np.array([q['mean'],q.ci_low,q.ci_high])*100
    coverage.errorbar(mean,3-j+offset,xerr=[[mean-lo],[hi-mean]],fmt=marker,ms=3,color=color,capsize=1.5,lw=.8)
  coverage.set(yticks=[3,2,1,0],yticklabels=[q[2] for q in configs],ylim=(-.5,4.2 if center=='NCEP' else 3.5),xlim=(-9,6),xticks=[-8,-4,0,4],xlabel='Coverage change (pp)')
  title(coverage,'cf'[i],'Coverage')
  if center=='NCEP':coverage.text(.98,.98,'Available = complete',ha='right',va='top',transform=coverage.transAxes,fontsize=8.5)
  change.text(.04,.98,'Same center by identity',ha='left',va='top',transform=change.transAxes,fontsize=8.5)
  for ax in (cloud,change,coverage):ax.tick_params(direction='out',pad=2,length=2.4,width=.65)
 handles=[Line2D([],[],marker='o',color=T,ms=3,lw=.8,label='MW'),Line2D([],[],marker='^',color=O,ms=3,lw=.8,label='CT'),Line2D([],[],marker='s',color='black',mfc='white',ms=4,ls='none',label='Mean'),Line2D([],[],marker='+',color='black',ms=5,ls='none',label='Truth'),Line2D([],[],color='#555555',lw=.8,label='Circle'),Line2D([],[],color='#555555',ls='--',lw=.8,label='Ellipse')]
 fig.legend(handles=handles,loc='center',bbox_to_anchor=(36/W,7/H),ncol=2,handlelength=1.1,columnspacing=.7,handletextpad=.4,labelspacing=.25)
 handles=[Line2D([],[],marker='s',color=T,ms=3.6,label='Center distance (d)'),Line2D([],[],marker='o',color=O,ms=3.6,label='Radius (r)')]
 fig.legend(handles=handles,loc='center',bbox_to_anchor=(98.5/W,7/H),ncol=1,handlelength=1.1,handletextpad=.4)
 handles=[Line2D([],[],marker='o',color=T,ms=3,label='Complete'),Line2D([],[],marker='D',color=P,ms=3,label='Available')]
 fig.legend(handles=handles,loc='center',bbox_to_anchor=(160.5/W,7/H),ncol=1,handlelength=1.1,handletextpad=.4)
 # Case provenance belongs to the caption; reserve the inter-row gap for axis labels.
 fig.canvas.draw()
 audit=require_matplotlib_panel_alignment(fig,axes=axes,panel_ids=list('abcdef'),row_groups=[['a','b','c'],['d','e','f']],column_groups=[['a','d'],['b','e'],['c','f']],json_out=out/'figure.alignment.json',overlay_svg=out/'figure.alignment.svg',tolerance_pt=1.5,require_panel_labels=True,strict=True,exemptions=[dict(panels=list('abcdef'),checks=['panel-width','horizontal-gutter'],reason='Explicit 44/37/35 mm columns and gaps reserve room for protocol labels; matching columns and row edges audited.')])
 for ext in ('pdf','svg','png'):fig.savefig(out/('Fig04_RangeDefinitions.'+ext),dpi=600)
 assert {p.name:SHA(p) for p in data.iterdir() if p.is_file()}==hashes
 (out/'render_complete.json').write_text(json.dumps(dict(status='PASS',input_sha256=hashes,source_sha256=SHA(__file__),size_mm=[W,H],font_floor_pt=8.5,alignment=audit['verdict'],numerical_changes=False,cloud_members=70,note='Original four-cloud display retained in SI.'),indent=2))
 plt.close(fig)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--font-dir',type=Path)
 a=p.parse_args();main(a.data,a.out,a.font_dir)
