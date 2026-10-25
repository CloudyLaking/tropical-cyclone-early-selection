"""Protected empirical evidence and nonduplicate SI panels; frozen data only."""
from pathlib import Path
import os,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from audit_panel_alignment import require_matplotlib_panel_alignment
from plot_fonts import configure_font
B=Path(__file__).resolve().parents[1];D=B/'data/evidence';O=Path(os.environ['REPRO_OUTPUT'])/'evidence';O.mkdir(parents=True,exist_ok=True)
T,S,K='#178F9F','#D48732','#252A30';R=['#B58C2B','#CB946A','#C8759F','#A382BC','#71549A'];markers=['o','v','s','^','D']
plt.rcParams.update({'font.family':configure_font(),'font.size':8,'axes.labelsize':8.5,'axes.titlesize':8.5,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.65,'axes.spines.top':True,'axes.spines.right':True,'legend.frameon':False,'figure.facecolor':'white'})
def choose(df,**kw):
 for k,v in kw.items():df=df[df[k]==v]
 return df

def title(ax,l,t):
 ax.text(0,1.06,l,transform=ax.transAxes,fontsize=10,weight='bold',va='bottom')
 ax.text(.06,1.06,t,transform=ax.transAxes,fontsize=8.5,va='bottom')

def save(fig,axes,name,cols):
 ids=list('abcdef')[:len(axes)];fig.canvas.draw()
 require_matplotlib_panel_alignment(fig,axes=axes,panel_ids=ids,row_groups=[ids[i:i+cols] for i in range(0,len(ids),cols) if len(ids[i:i+cols])>1],column_groups=[ids[i::cols] for i in range(cols)] if len(ids)>cols else [],json_out=O/(name+'.alignment.json'),strict=True,require_panel_labels=True)
 for ext in ['pdf','svg','png']:fig.savefig(O/(name+'.'+ext),dpi=600)
 plt.close(fig)

h=pd.read_csv(B/'data/si_readability/member_histories.csv.gz');r=pd.read_csv(D/'rank_806_summary.csv')
f=plt.figure(figsize=(180/25.4,134/25.4));axes=[f.add_axes([x/180,y/134,65/180,37/134]) for y in [79,24] for x in [16,109]]
counts={}
for row,center in enumerate(['ECMWF','NCEP']):
 a,b=axes[2*row:2*row+2];date=['2015-05-12T12:00:00','2015-10-18T00:00:00'][row]
 z=choose(h,center=center,base=date,scoring='raw');counts[center]=z.member.nunique()
 for mid,q in z.groupby('member'):
  q=q.sort_values('lead');g=min(4,int(q.early_fraction.iloc[0]*5));a.plot(q.lead,100*q.rank_fraction,color=R[g],ls='-' if q.selected.iloc[0] else '--',lw=.55,alpha=.9)
 a.set(xlim=(24,120),ylim=(100,0),xticks=[24,48,72,96,120],yticks=[0,25,50,75,100],xlabel='Lead time (h)',ylabel='Raw member-error rank (%)');title(a,'ac'[row],center+' · '+date[:10]+' '+date[11:13]+' UTC')
 z=choose(r,center=center,lead=120).set_index('metric');assert set(z.cases)=={300 if row==0 else 506}
 b.axhline(20,color='#8B949C',ls=':',lw=.8)
 for k in range(5):
  q=z.loc[[f'm{k}{j}' for j in range(5)]];y=q['mean'].to_numpy()*100
  b.errorbar(np.arange(5),y,yerr=[y-q.ci_low*100,q.ci_high*100-y],color=R[k],marker=markers[k],ms=3,lw=.9,capsize=1.5)
 b.set(xlim=(-.3,4.3),ylim=(10,34),xticks=range(5),xticklabels=['Best','2','3','4','Worst'],xlabel='Aligned 120-h error group',ylabel='Later-group members (%)');title(b,'bd'[row],center+' · matched transitions')
 b.text(.03,.97,f'{300 if row==0 else 506} initializations / {64 if row==0 else 100} storms',transform=b.transAxes,fontsize=7.5,va='top')
f.text(.5,130/134,'Raw 24-h groups; left: raw rank histories · right: aligned 120-h verification',ha='center',fontsize=8)
f.legend(handles=[Line2D([],[],color=c,marker=m,label=l,lw=.8,ms=3) for c,m,l in zip(R,markers,['Early best','2','3','4','Early worst'])],loc='center',bbox_to_anchor=(.5,10/134),ncol=5,columnspacing=1.2)
f.text(.5,2/134,'Histories: solid = selected; dashed = excluded. 20%: no association between equal-mass quintiles.',ha='center',fontsize=7.5)
save(f,axes,'evidence_rank_histories_and_transitions',2)
(O/'rank_receipt.json').write_text(json.dumps({'member_counts':counts,'transition_cases':[300,506],'new_statistics':False},indent=2))

examples=[]
for center,name in [('ECMWF','illustrative_track.npz'),('NCEP','illustrative_track_ncep.npz')]:
 z=np.load(D/name,allow_pickle=False);truth=z['truth_xy'];tracks=z['p']+truth[:,None,:];origin=truth[2];examples.append((center,truth-origin,tracks-origin,z['selected']))
allxy=np.concatenate([q.reshape(-1,2) for _,tr,p,ix in examples for q in [tr[2:],p[2:]]]);lo=allxy.min(0);hi=allxy.max(0);cx=(hi+lo)/2;span=max(hi-lo)*1.12
f=plt.figure(figsize=(180/25.4,105/25.4));axes=[f.add_axes([x/180,23/105,65/180,65/105]) for x in [17,110]]
for a,(center,tr,p,ix),letter in zip(axes,examples,'ab'):
 for j in range(p.shape[1]):a.plot(p[2:,j,0],p[2:,j,1],color=S if j in ix else '#AEB7BC',alpha=.38,lw=.45)
 for arr,c,m,label in [(p.mean(1),T,'s','Full mean'),(p[:,ix].mean(1),S,'D','Early half')]:a.plot(arr[2:,0],arr[2:,1],color=c,lw=1.6,marker=m,markevery=[0,4,8],ms=4,label=label)
 a.plot(tr[2:,0],tr[2:,1],color=K,lw=1.5)
 for k,hour in [(2,24),(6,72),(10,120)]:
  a.plot(*tr[k],marker='+',ms=6,color=K);a.annotate(str(hour)+' h',tr[k],xytext=(18,-7) if hour==120 else (5,-7 if hour==24 else 4),textcoords='offset points',fontsize=8)
 a.set(xlim=(cx[0]-span/2,cx[0]+span/2),ylim=(cx[1]-span/2,cx[1]+span/2),xlabel='East of observed 24-h position (km)',ylabel='North of observed 24-h position (km)');a.set_aspect('equal');title(a,letter,f'{center} · {p.shape[1]} members / {len(ix)} selected')
f.text(.5,100/105,'Hagibis (2019) · initialized 7 October 2019, 12 UTC',ha='center',fontsize=9)
f.legend(handles=[Line2D([],[],color=T,marker='s',label='Full mean'),Line2D([],[],color=S,marker='D',label='Early half'),Line2D([],[],color=K,marker='+',label='Observed')],loc='center',bbox_to_anchor=(.5,5/105),ncol=3)
save(f,axes,'evidence_hagibis_aligned_tracks',2)

# Old S3c,f: sample sensitivity only.
p=pd.read_csv(D/'performance_summary.csv');f=plt.figure(figsize=(180/25.4,69/25.4));axes=[f.add_axes([x/180,20/69,65/180,38/69]) for x in [16,109]]
for a,center,l in zip(axes,['ECMWF','NCEP'],'ab'):
 for cohort,c,m,ls in [('available','#8B949C','o',':'),('common144','#8062A6','^','--'),('common168',S,'D','-')]:
  q=choose(p,center=center,cohort=cohort,metric='selection_after_translation_km').sort_values('lead')
  a.errorbar(q.lead,q.value,yerr=[q.value-q.lo,q.hi-q.value],color=c,marker=m,ls=ls,ms=3,capsize=1.5,lw=1,label={'available':'Available at lead','common144':'Fixed 144 h','common168':'Fixed 168 h'}[cohort])
 a.axhline(0,color=K,lw=.6);a.set(xlabel='Lead time (h)',ylabel='Additional selection gain (km)',xticks=[48,72,96,120,144,168]);title(a,l,center)
f.legend(*axes[0].get_legend_handles_labels(),loc='center',bbox_to_anchor=(.5,4/69),ncol=3)
save(f,axes,'S_sample_sensitivity',2)

# Old S18a-d: complete-cohort leadwise circle and ellipse coverage.
d=pd.read_csv(B/'data/protocol/results/summary.csv');f=plt.figure(figsize=(180/25.4,121/25.4));axes=[f.add_axes([x/180,y/121,65/180,37/121]) for y in [74,20] for x in [16,109]]
for row,(shape,metric) in enumerate([('Circle','in70'),('Ellipse','ellipse_in70')]):
 for col,center in enumerate(['ECMWF','NCEP']):
  a=axes[row*2+col]
  for group,c,m in [('random','#8B949C','o'),('selected',S,'s')]:
   for protocol,ls,filled in [('A','-',True),('B','--',False)]:
    q=choose(d,cohort='main806',center=center,comparison=protocol+'_'+group,metric=metric).sort_values('lead')
    a.fill_between(q.lead,100*q.ci_low,100*q.ci_high,color=c,alpha=.07,lw=0);a.plot(q.lead,100*q['mean'],color=c,marker=m,ls=ls,mfc=c if filled else 'white',ms=3,lw=1.1,label=group.capitalize()+' · '+('MW' if filled else 'CT'))
  a.set(xlim=(44,124),ylim=(20,100),xticks=[48,72,96,120],xlabel='Lead time (h)',ylabel='Truth covered (%)');title(a,'abcd'[row*2+col],center+' · '+shape)
f.legend(*axes[0].get_legend_handles_labels(),loc='center',bbox_to_anchor=(.5,3/121),ncol=4)
save(f,axes,'S_leadwise_coverage',2)

# Complete tails of the retained-fraction distribution; no whisker-window crop.
st=pd.read_csv(D/'fraction_storms.csv');st=st[st.lead==120];f=plt.figure(figsize=(180/25.4,105/25.4));axes=[f.add_axes([17/180,y/105,158/180,29/105]) for y in [63,20]]
for a,metric,l in zip(axes,['ATE','CTE'],'ab'):
 for center,offset,c in [('ECMWF',-2,'#5675A7'),('NCEP',2,'#B46B87')]:
  for frac in range(10,101,10):
   q=choose(st,center=center,fraction=frac).sort_values('sid');assert len(q)==(64 if center=='ECMWF' else 100)
   vals=q[metric].to_numpy();x=frac+offset
   a.boxplot(vals,positions=[x],widths=3,whis=(5,95),showfliers=True,patch_artist=True,manage_ticks=False,boxprops={'facecolor':c+'35','edgecolor':c,'lw':.7},medianprops={'color':c,'lw':1},whiskerprops={'color':c,'lw':.65},capprops={'color':c,'lw':.65},flierprops={'marker':'.','markersize':2,'markeredgewidth':0,'markerfacecolor':c,'alpha':.5})
   a.plot(x,vals.mean(),'+',color=K,ms=3.2,mew=.8)
 a.axhline(0,color='#777777',lw=.65);a.set(xlim=(4,106),xticks=range(10,101,10),ylabel=metric+' (km)');title(a,l,'120-h '+metric+' across retained fractions')
axes[1].set_xlabel('Requested retention (%)')
f.legend(handles=[Line2D([],[],color=c,lw=3,label=l) for c,l in [('#5675A7','ECMWF (64 storms)'),('#B46B87','NCEP (100 storms)')]]+[Line2D([],[],color=K,marker='+',ls='none',label='Storm-equal mean')],loc='center',bbox_to_anchor=(.5,4/105),ncol=3)
save(f,axes,'S17_retention_fullrange',1)
