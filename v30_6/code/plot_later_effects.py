"""Display frozen forecast effects without new statistics."""
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
T,O='#0095A6','#E98B2A'
def main(a):
 a.out.mkdir(parents=True,exist_ok=False)
 FONT=configure_font(a.fonts)
 plt.rcParams.update({'font.family':FONT,'font.size':8.5,'axes.labelsize':8.5,'axes.titlesize':8.5,'xtick.labelsize':8.5,'ytick.labelsize':8.5,'legend.fontsize':8.5,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.65,'axes.spines.top':True,'axes.spines.right':True,'legend.frameon':False})
 perf=pd.read_csv(a.data/'performance_dense.csv');storms=pd.read_csv(a.data/'fixed_storms.csv');summary=pd.read_csv(a.data/'fixed_summary.csv')
 width_mm=183;W,H=width_mm,114;fig=plt.figure(figsize=(W/25.4,H/25.4))
 axes=[fig.add_axes([x/W,y/H,43/W,29/H]) for y in [72,22] for x in [15,75,135]]
 fig.text(66.5/W,111/H,'Fixed 168-h cohort · great-circle distance',ha='center',va='top',fontsize=8.5)
 fig.text(156.5/W,111/H,'Complete 120 h · projected errors',ha='center',va='top',fontsize=8.5)
 used=[]; paired_labels=[]
 geometry=pd.read_csv(a.data/'geometry_cases.csv')
 def title(ax,l,s):
  ax.text(0,1.08,l,transform=ax.transAxes,fontweight='bold',fontsize=9,va='bottom')
  ax.text(.09,1.08,s,transform=ax.transAxes,fontsize=8.5,va='bottom')
 for i,center in enumerate(['ECMWF','NCEP']):
  error,gain,dist=axes[3*i:3*i+3]
  p=perf[(perf.center==center)&(perf.cohort=='common168')&(perf.lead%24==0)]
  assert set(p.n_cases)==({100} if i==0 else {195})
  for method,color,marker,ls,filled in [('raw_full',T,'s','--',False),('raw_selected',O,'D','--',False),('translated_full',T,'s','-',True),('translated_selected',O,'D','-',True)]:
   q=p[p.metric==method].sort_values('lead');assert q.lead.tolist()==list(range(0,169,24))
   if method!='raw_full':q=q[q.lead>=24]
   error.plot(q.lead,q.value,color=color,marker=marker,ls=ls,lw=1,ms=3.7,mfc=color if filled else 'white',mew=.8,zorder=4 if filled else 2)
   used.extend(q.to_dict('records'))
  title(error,'ad'[i],center+' · error')
  error.set(xlim=(-4,174),ylim=(0,750),xticks=[0,24,72,120,168],yticks=[0,250,500,750],xlabel='Forecast lead (h)',ylabel='Position error (km)')
  error.text(.21,.96,'24-h update',transform=error.transAxes,va='top',fontsize=8.5)
  error.text(.21,.77,('100 initializations\n21 storms' if i==0 else '195 initializations\n41 storms'),transform=error.transAxes,va='top',fontsize=8.5)
  for ax in [error,gain]:
   ax.axvspan(0,24,color='#D9DEE1',alpha=.35,zorder=0)
   ax.axvline(24,color='#888888',ls=':',lw=.7,zorder=1)
  gain.axhline(0,color='#777777',lw=.65)
  for metric,color,marker in [('translation_full_km',T,'s'),('selection_after_translation_km',O,'D')]:
   q=p[(p.metric==metric)&(p.lead>=24)].sort_values('lead');assert len(q)==7
   post=q[q.lead>24];anchor=q[q.lead==24]
   gain.errorbar(post.lead,post.value,yerr=[post.value-post.lo,post.hi-post.value],color=color,marker=marker,lw=1,ms=3.7,capsize=2,elinewidth=.65)
   gain.errorbar(anchor.lead,anchor.value,yerr=[(anchor.value-anchor.lo).to_numpy(),(anchor.hi-anchor.value).to_numpy()],color='#999999',marker=marker,mfc='white',ls='none',lw=.65,ms=3.7,capsize=2,elinewidth=.65)
   used.extend(q.to_dict('records'))
  gain.text(.24,.97,'24 h: imposed',transform=gain.transAxes,color='#666666',va='top',fontsize=8.5)
  title(gain,'be'[i],center+' · improvement')
  gain.set(xlim=(-4,174),ylim=(-56,80),xticks=[0,24,72,120,168],yticks=[-40,0,40,80],xlabel='Forecast lead (h)',ylabel='Error reduction (km)')
  title(dist,'cf'[i],center+' · 120 h')
  dist.axhline(0,color='#888888',lw=.65)
  for k,metric in enumerate(['ATE','CTE']):
   for j,method in enumerate(['full','early']):
    q=storms[(storms.center==center)&(storms.lead==120)&(storms.method==method)].sort_values('sid')
    assert len(q)==(64 if i==0 else 100) and q.sid.is_unique
    vals=q[metric].to_numpy();assert vals.min()>-1100 and vals.max()<1400
    color=[T,O][j];x=k+[-.16,.16][j]
    jitter=np.array([(int(hashlib.sha256(str(s).encode()).hexdigest()[:8],16)%1000/999-.5)*.2 for s in q.sid])
    dist.scatter(x+jitter,vals,s=4.5,color=color,alpha=.3,linewidths=0)
    dist.boxplot(vals,positions=[x],widths=.25,whis=(5,95),showfliers=False,patch_artist=True,manage_ticks=False,boxprops={'facecolor':color+'25','edgecolor':color,'lw':.8},medianprops={'color':color,'lw':1},whiskerprops={'color':color,'lw':.8},capprops={'color':color,'lw':.8})
    mean=float(vals.mean());ref=summary[(summary.center==center)&(summary.lead==120)&(summary.method==method)&(summary.metric==metric)].iloc[0]['mean'];assert abs(mean-ref)<1e-7
    dist.plot(x,mean,'+',color='#252525',ms=5,mew=1)
  dist.set(xticks=[0,1],xticklabels=['ATE','CTE'],xlim=(-.55,1.55),ylim=(-1100,1400),yticks=[-1000,0,1000],ylabel='Storm-mean error (km)')
  dist.text(.04,.96,('300 initializations\n64 storms' if i==0 else '506 initializations\n100 storms'),transform=dist.transAxes,va='top')
  contrast=summary[(summary.center==center)&(summary.lead==120)&(summary.method=='early_minus_full')&(summary.metric=='distance')].iloc[0]
  mean_gain=-float(contrast['mean']);lo=-float(contrast.ci_high);hi=-float(contrast.ci_low)
  geo=geometry[(geometry.center==center)&(geometry.lead==120)]
  assert len(geo)==(300 if i==0 else 506) and abs(geo.weight.sum()-1)<1e-10
  improve=100*float((geo.improved*geo.weight).sum())
  assert round(improve,2)==(55.40 if i==0 else 53.68)
  assert abs(float((geo.gain*geo.weight).sum())-mean_gain)<1e-7
  label=f'{improve:.1f}% improve; mean distance gain\n{mean_gain:.1f} [{lo:.1f}, {hi:.1f}] km (95% CI)'.replace('-', '−')
  dist.set_xlabel(label,labelpad=4,linespacing=1.3)
  paired_labels.append(dict(center=center,lead=120,improved_percent=improve,mean_gain=mean_gain,ci_low=lo,ci_high=hi,baseline='aligned full minus aligned early half',cohort_cases=len(geo)))
  for ax in [error,gain,dist]:ax.tick_params(length=2.4,width=.65,pad=2)
 # Figure-level legends are tied to column roles, not inter-row dead space.
 handles=[Line2D([],[],color=T,marker='s',ms=3.5,label='Full / alignment gain'),Line2D([],[],color=O,marker='D',ms=3.5,label='Selected / added-selection gain'),Line2D([],[],color='#555555',ls='--',marker='o',mfc='white',ms=3.5,label='Raw'),Line2D([],[],color='#555555',ls='-',marker='o',ms=3.5,label='Aligned')]
 fig.legend(handles=handles,loc='center',bbox_to_anchor=(.51,5/H),ncol=4,columnspacing=.9,handletextpad=.4,handlelength=1.4)
 fig.canvas.draw()
 alignment=require_matplotlib_panel_alignment(fig,axes=axes,panel_ids=list('abcdef'),row_groups=[list('abc'),list('def')],column_groups=[list('ad'),list('be'),list('cf')],json_out=a.out/'alignment.json',overlay_svg=a.out/'alignment.svg',tolerance_pt=1.5,strict=True,require_panel_labels=True)
 fig.savefig(a.out/'Fig02_EmpiricalResponse.pdf',dpi=600)
 fig.savefig(a.out/'Fig02_EmpiricalResponse.svg',dpi=600)
 fig.savefig(a.out/'Fig02_EmpiricalResponse.png',dpi=600)
 pd.DataFrame(paired_labels).to_csv(a.out/'paired_labels.csv',index=False)
 pd.DataFrame(used).to_csv(a.out/'plotted_performance_rows.csv',index=False)
 hashes={n:hashlib.sha256((a.data/n).read_bytes()).hexdigest() for n in ['performance_dense.csv','fixed_storms.csv','fixed_summary.csv','geometry_cases.csv']}
 (a.out/'acceptance.json').write_text(json.dumps({'status':'PASS','input_sha256':hashes,'size_mm':[W,H],'font_floor_pt':8.5,'alignment':alignment['verdict'],'new_forecast_methods':False,'added_verification_leads':[],'display_revision':'Paired full-cohort gain and storm-weighted improvement proportion; unchanged frozen intervals','plotted_grid_hours':list(range(0,169,24)),'fixed168_sample':{'ECMWF':[100,21],'NCEP':[195,41]},'120h_distribution_sample':{'ECMWF':[300,64],'NCEP':[506,100]}},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--fonts',type=Path);main(p.parse_args())
