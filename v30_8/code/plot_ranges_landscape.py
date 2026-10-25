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
from matplotlib.ticker import MaxNLocator
from matplotlib.cbook import boxplot_stats
from audit_panel_alignment import require_matplotlib_panel_alignment
SHA=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
T,O,P='#178F9F','#D48732','#8062A6'
def main(data,out,fonts):
 out.mkdir(parents=True,exist_ok=False)
 hashes={p.name:SHA(p) for p in data.iterdir() if p.is_file()}
 FONT=configure_font(fonts)
 plt.rcParams.update({'font.family':FONT,'font.size':8.5,'axes.labelsize':8.5,'axes.titlesize':8.5,'xtick.labelsize':8.5,'ytick.labelsize':8.5,'legend.fontsize':8.5,'axes.linewidth':.65,'pdf.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white','axes.facecolor':'white','axes.spines.top':True,'axes.spines.right':True,'legend.frameon':False})
 read=lambda n:pd.read_csv(data/n)
 m=read('source_cloud_members.csv');b=read('source_region_boundaries.csv');g=read('source_region_geometry.csv');dr=read('source_center_radius_changes.csv');cov=read('source_paired_coverage.csv')
 assert len(m)==70 and len(cov)==16 and len(dr)==8
 assert json.loads((data/'margin_acceptance.json').read_text())['status']=='PASS'
 # Distribution units are storms. Reuse saved case scores; do not rescore tracks.
 source_cases=data.parent/'protocol/results/cases.csv.gz';source_case_hash=SHA(source_cases)
 case_scores=pd.read_csv(source_cases,float_precision='round_trip');case_scores=case_scores[case_scores.lead==120]
 storm_values={};distribution_rows=[]
 for cohort in ['main806','available']:
  subset=case_scores[case_scores.main806] if cohort=='main806' else case_scores
  for center in ['ECMWF','NCEP']:
   for protocol in ['A','B']:
    part=subset[(subset.center==center)&(subset.protocol==protocol)&subset['group'].isin(['selected','random'])]
    for metric in ['distance','r70','in70','ellipse_in70']:
     paired=part.pivot(index=['sid','code'],columns='group',values=metric)
     assert not paired[['selected','random']].isna().any().any()
     vals=(paired['selected']-paired['random']).groupby(level='sid').mean()
     assert len(vals)==({'ECMWF':64,'NCEP':100}[center] if cohort=='main806' else {'ECMWF':69,'NCEP':100}[center])
     storm_values[cohort,center,protocol,metric]=vals.to_numpy()
     for sid,value in vals.items():distribution_rows.append(dict(cohort=cohort,center=center,protocol=protocol,metric=metric,sid=sid,value=value))
 radius_pool=np.concatenate([v for (cohort,center,protocol,metric),v in storm_values.items() if cohort=='main806' and metric in ['distance','r70']])
 span=np.ptp(radius_pool);full_distance_limits=(radius_pool.min()-.04*span,radius_pool.max()+.04*span)
 def visible_limits(metrics,summary,scale=1):
  bounds=[0.]
  for (cohort,center,protocol,metric),values in storm_values.items():
   if metric not in metrics or (metric in ['distance','r70'] and cohort!='main806'):continue
   stat=boxplot_stats(values*scale,whis=(5,95))[0]
   bounds.extend([stat['whislo'],stat['whishi']])
  bounds.extend((summary[['ci_low','ci_high']].to_numpy()*scale).ravel())
  low,high=min(bounds),max(bounds);pad=.07*(high-low)
  step=10 if scale==1 else 5
  return (np.floor((low-pad)/step)*step,np.ceil((high+pad)/step)*step)
 distance_limits=visible_limits(['distance','r70'],dr)
 coverage_limits=visible_limits(['in70','ellipse_in70'],cov,100)
 display_rows=[];edge_labels=[]
 def mark_outside(ax,arrays,limits,panel):
  values=np.concatenate(arrays);left=int((values<limits[0]).sum());right=int((values>limits[1]).sum())
  # Totals count storm entries across the separate boxes in this panel.
  edge_labels.append(ax.text(.03,.97,f'← {left}',ha='left',va='top',transform=ax.transAxes,fontsize=8,color='#555555'))
  edge_labels.append(ax.text(.97,.97,f'{right} →',ha='right',va='top',transform=ax.transAxes,fontsize=8,color='#555555'))
  display_rows.append(dict(panel=panel,limits=list(limits),outside_left=left,outside_right=right,entries=len(values)))
 flier_restore=[]
 def box(ax,values,y,color,width=.20,available=False):
  stats=boxplot_stats(values,whis=(5,95))[0];all_fliers=stats['fliers'].copy()
  limits=distance_limits if width==.20 else coverage_limits
  stats['fliers']=all_fliers[(all_fliers>=limits[0])&(all_fliers<=limits[1])]
  artists=ax.bxp([stats],positions=[y],vert=False,widths=width,showfliers=True,patch_artist=True,manage_ticks=False,
   boxprops={'facecolor':'white' if available else color+'35','edgecolor':color,'lw':.7,'linestyle':'--' if available else '-'},
   medianprops={'color':color,'lw':1},whiskerprops={'color':color,'lw':.6},capprops={'color':color,'lw':.6},
   flierprops={'marker':'.','markersize':1.8,'markeredgewidth':0,'markerfacecolor':color,'alpha':.35})
  flier_restore.append((artists['fliers'][0],all_fliers,y))

 W,H=180,124
 fig=plt.figure(figsize=(W/25.4,H/25.4))
 axes=[fig.add_axes([x/W,y/H,43/W,30/H]) for y in (76,27) for x in (14,74,134)]
 def title(ax,l,t):
  ax.text(0,1.07,l,transform=ax.transAxes,fontsize=10,fontweight='bold',va='bottom')
  ax.text(.10,1.07,t,transform=ax.transAxes,fontsize=9.5,va='bottom')
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
  cloud.set(xlim=(-lim*43/30,lim*43/30),ylim=(-lim,lim),xticks=[-400,0,400],yticks=[-400,0,400],xlabel='East of selected mean (km)',ylabel='North of selected mean (km)')
  cloud.set_aspect('equal',adjustable='box')
  cloud.text(.03,.97,'Circle: covered if d ≤ r',transform=cloud.transAxes,ha='left',va='top',fontsize=8.5)
  title(cloud,'ad'[i],center+' · same mean')
  change_arrays=[];coverage_arrays=[]
  change.axvline(0,ymax=.82,color='#777777',ls='--',lw=.65)
  for j,protocol in enumerate(('A','B')):
   for metric,color,offset,marker in [('distance',T,.12,'s'),('r70',O,-.12,'D')]:
    q=dr[(dr.center==center)&(dr.protocol==protocol)&(dr.metric==metric)].iloc[0]
    values=storm_values['main806',center,protocol,metric];assert abs(values.mean()-q['mean'])<1e-8
    change_arrays.append(values)
    box(change,values,1-j+offset,color,width=.20)
    change.errorbar(q['mean'],1-j+offset,xerr=[[q['mean']-q.ci_low],[q.ci_high-q['mean']]],fmt=marker,ms=3.6,color='#252A30',mfc=color,capsize=1.8,lw=1.2,zorder=6)
  change.set(yticks=[1,0],yticklabels=['MW','CT'],ylim=(-.5,1.5),xlim=distance_limits,xlabel='Selected − same-size\nrandom (km)')
  change.xaxis.set_major_locator(MaxNLocator(3))
  title(change,'be'[i],'Center / radius')
  coverage.axvline(0,ymax=.82 if center=='NCEP' else 1,color='#777777',ls='--',lw=.65)
  configs=[('A','in70','MW circle'),('B','in70','CT circle'),('A','ellipse_in70','MW ellipse'),('B','ellipse_in70','CT ellipse')]
  cohorts=[('main806','#33373C',.18,'o'),('available','#33373C',-.18,'D')] if center=='ECMWF' else [('main806','#33373C',0,'o')]
  for cohort,color,offset,marker in cohorts:
   for j,(protocol,metric,label) in enumerate(configs):
    q=cov[(cov.center==center)&(cov.cohort==cohort)&(cov.comparison==protocol+'_selected_minus_random')&(cov.metric==metric)].iloc[0]
    mean,lo,hi=np.array([q['mean'],q.ci_low,q.ci_high])*100
    color=T if protocol=='A' else O
    values=100*storm_values[cohort,center,protocol,metric];assert abs(values.mean()-mean)<1e-8
    coverage_arrays.append(values)
    box(coverage,values,3-j+offset,color,width=.24,available=(cohort=='available'))
    coverage.errorbar(mean,3-j+offset,xerr=[[mean-lo],[hi-mean]],fmt=marker,ms=3,color='#252A30',mfc='white' if cohort=='available' else color,capsize=1.8,lw=1.2,zorder=6)
  coverage.set(yticks=[3,2,1,0],yticklabels=[q[2] for q in configs],ylim=(-.5,3.85),xlim=coverage_limits,xlabel='Selected − random coverage\n(percentage points)')
  title(coverage,'cf'[i],'Coverage')
  coverage.xaxis.set_major_locator(MaxNLocator(3))
  mark_outside(change,change_arrays,distance_limits,'be'[i])
  mark_outside(coverage,coverage_arrays,coverage_limits,'cf'[i])
  for ax in (cloud,change,coverage):ax.tick_params(direction='out',pad=2,length=2.4,width=.65)
 handles=[Line2D([],[],marker='o',color=T,ms=3,lw=.8,label='MW'),Line2D([],[],marker='^',color=O,ms=3,lw=.8,label='CT'),Line2D([],[],marker='s',color='black',mfc='white',ms=4,ls='none',label='Mean'),Line2D([],[],marker='+',color='black',ms=5,ls='none',label='Truth'),Line2D([],[],color='#555555',lw=.8,label='Circle'),Line2D([],[],color='#555555',ls='--',lw=.8,label='Ellipse')]
 fig.legend(handles=handles,loc='center',bbox_to_anchor=(35.5/W,11/H),ncol=2,handlelength=1.1,columnspacing=.7,handletextpad=.4,labelspacing=.25)
 handles=[Line2D([],[],marker='s',color=T,ms=3.6,label='Center error d'),Line2D([],[],marker='D',color=O,ms=3.6,label='Member radius r')]
 fig.legend(handles=handles,loc='center',bbox_to_anchor=(95.5/W,11/H),ncol=1,handlelength=1.1,handletextpad=.4)
 handles=[Line2D([],[],marker='o',color='#33373C',ms=3,label='Complete'),Line2D([],[],marker='D',color='#33373C',mfc='white',ms=3,label='Available')]
 fig.legend(handles=handles,loc='center',bbox_to_anchor=(155.5/W,11/H),ncol=1,handlelength=1.1,handletextpad=.4)
 pd.DataFrame(distribution_rows).to_csv(out/'plotted_storm_changes.csv',index=False)
 assert SHA(source_cases)==source_case_hash
 # Case provenance belongs to the caption; reserve the inter-row gap for axis labels.
 for x,s in [(35.5,'Same mean, different clouds'),(95.5,'Center and radius changes'),(155.5,'Coverage changes')]:fig.text(x/W,118/H,s,ha='center',va='top',fontsize=8.5)
 fig.text(.5,2.7/H,'Boxes: storms · Black bars: mean 95% CI · Edge counts: off-scale entries',ha='center',va='center',fontsize=8)
 fig.canvas.draw()
 audit=require_matplotlib_panel_alignment(fig,axes=axes,panel_ids=list('abcdef'),row_groups=[['a','b','c'],['d','e','f']],column_groups=[['a','d'],['b','e'],['c','f']],json_out=out/'figure.alignment.json',overlay_svg=out/'figure.alignment.svg',tolerance_pt=1.5,require_panel_labels=True,strict=True)
 for ext in ('pdf','svg','png'):fig.savefig(out/('Fig04_RangeDefinitions.'+ext),dpi=600)
 (out/'display_limits.json').write_text(json.dumps(dict(rule='Common limits across centers enclose all 5th/95th-percentile whiskers and original mean confidence intervals, plus 7% padding rounded outward.',panels=display_rows),indent=2))
 for label in edge_labels:label.set_visible(False)
 for artist,values,y in flier_restore:artist.set_data(values,np.full(len(values),y))
 for ax in [axes[1],axes[4]]:ax.set_xlim(full_distance_limits)
 for ax in [axes[2],axes[5]]:ax.set_xlim(-105,105)
 fig.texts[-1].set_text('Full range: all storm entries · Black bars: mean 95% CI')
 for ext in ('pdf','svg','png'):fig.savefig(out/('Fig04_RangeDefinitions_FullRange.'+ext),dpi=600)
 assert {p.name:SHA(p) for p in data.iterdir() if p.is_file()}==hashes
 (out/'render_complete.json').write_text(json.dumps(dict(status='PASS',input_sha256=hashes,source_sha256=SHA(__file__),size_mm=[W,H],font_floor_pt=8.5,alignment=audit['verdict'],numerical_changes=False,distribution_unit='storm-level selected-minus-random changes from frozen case scores',box_whiskers='5th and 95th percentiles; off-scale entry counts shown; all fliers retained in full-range companion',source_case_sha256=SHA(source_cases),cloud_members=70,note='Original four-cloud display retained in SI.'),indent=2))
 plt.close(fig)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--font-dir',type=Path)
 a=p.parse_args();main(a.data,a.out,a.font_dir)
