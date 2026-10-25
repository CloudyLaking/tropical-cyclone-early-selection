from pathlib import Path
import json,os
from plot_fonts import configure_font
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt, font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Arc
from audit_panel_alignment import require_matplotlib_panel_alignment

B=Path(__file__).resolve().parent.parent; D=B/'data/evidence'; O=Path(os.environ['REPRO_OUTPUT'])/'evidence'; O.mkdir(exist_ok=True,parents=True)
FONT=configure_font()
plt.rcParams.update({'font.family':FONT,'font.size':8,'axes.labelsize':8,'axes.titlesize':8.5,
    'xtick.labelsize':7.5,'ytick.labelsize':7.5,'axes.linewidth':.65,'axes.edgecolor':'#333333',
    'pdf.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white','axes.facecolor':'white',
    'legend.frameon':False,'lines.linewidth':1.1,'savefig.facecolor':'white','mathtext.fontset':'stix'})
C={'full':'#178F9F','early':'#D48732','selected':'#D48732','path':'#9457BA','terminal':'#D0A018','random':'#777C86'}
M={'full':'s','early':'D','selected':'D','path':'^','terminal':'v','random':'o'}
def select(df,**kw):
    m=np.ones(len(df),bool)
    for k,v in kw.items():m&=df[k]==v
    return df[m]
def row(df,**kw):
    x=select(df,**kw);assert len(x)==1,(kw,len(x));return x.iloc[0]
def grid(h=132,bottom=19):
    fig=plt.figure(figsize=(7.2047244094,h/25.4)); axes=[]
    for r in range(2):
        for c in range(3):axes.append(fig.add_axes([(16+c*57)/183,(bottom+(1-r)*62)/h,44/183,44/h]))
    return fig,axes
def title(ax,letter,t):
    bounds=ax.get_position();fw,fh=ax.figure.get_size_inches()*25.4
    y=1+2/(bounds.height*fh)
    ax.text(0,y,letter,transform=ax.transAxes,weight='bold',fontsize=9,va='bottom')
    ax.text(4.2/(bounds.width*fw),y,t,transform=ax.transAxes,fontsize=8.2,va='bottom')
def handle(k,label,ls='-'):
    return Line2D([],[],color=C[k],marker=M[k],ms=3.4,lw=1.1,ls=ls,label=label)
def save(fig,axs,name):
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig,json_out=O/(name+'_alignment.json'),axes=axs,panel_ids=list('abc'),
      row_groups=[list('abc')],column_groups=[],
      tolerance_pt=1.5,gutter_tolerance_pt=1.5,require_panel_labels=True,strict=True)
    for ext in ['pdf','svg','png']:fig.savefig(O/(name+'.'+ext),dpi=600)
    plt.close(fig)

parts=pd.read_csv(D/'lag_decomposition_TableS3.csv')
fig=plt.figure(figsize=(180/25.4,79/25.4))
positions=[(14,19,44,44),(73,19,44,44),(132,19,44,44)]
axs=[fig.add_axes([x/180,y/79,w/180,h/79]) for x,y,w,h in positions]
fig.text(.5,76/79,'Displacement geometry and mean-lag contributions',ha='center',va='top',fontsize=9)
ax=axs[0];title(ax,'a','Hagibis · 108–120 h')
z=np.load(D/'illustrative_track.npz'); truth=z['truth_xy']; pos=z['p']+truth[:,None,:]
f=pos[-1]-pos[-2]; obs=truth[-1]-truth[-2]; lo=np.linalg.norm(obs); u=obs/lo; v=np.array([-u[1],u[0]])
rot=np.stack([u,v],axis=1); fr=f@rot; fm=fr.mean(axis=0); sm=fr[z['selected']].mean(axis=0)
L=np.linalg.norm(fr,axis=1).mean(); s=np.linalg.norm(fm)
terms=np.array([[L-lo,0],[-(L-s),0],fm-np.array([s,0])])
assert np.allclose(terms.sum(axis=0),fm-np.array([lo,0]),atol=1e-10)
# Every actual member displacement is shown; these are vectors, not positions.
for vec in fr:ax.plot([0,vec[0]],[0,vec[1]],color='#BEC5CA',alpha=.6,lw=.55,zorder=1)
ax.scatter(fr[:,0],fr[:,1],s=4,color='#98A4AA',alpha=.5,lw=0,zorder=2)
for end,col,mark in [(fm,C['full'],'s'),(sm,C['early'],'D'),(np.array([lo,0]),'#222222','+')]:
    ax.annotate('',xy=end,xytext=(0,0),arrowprops={'arrowstyle':'->','color':col,'lw':1.4,'shrinkA':0,'shrinkB':0},zorder=4)
    ax.plot(*end,marker=mark,ms=3.6,color=col,zorder=5)
for rad,style in [(L,'--'),(s,':')]:ax.add_patch(Arc((0,0),2*rad,2*rad,theta1=-18,theta2=50,color='#62686C',lw=.8,ls=style,zorder=3))
ax.text(.04,.96,f'ECMWF · 51 members',transform=ax.transAxes,va='top',fontsize=7.5)
ax.text(.04,.24,f'L = {L:.0f} km\ns = {s:.0f} km\nobs = {lo:.0f} km',transform=ax.transAxes,va='top',fontsize=7.5,linespacing=1.15)
xmin=min(0,fr[:,0].min()); xmax=max(lo,fr[:,0].max())
ymin=min(0,fr[:,1].min()); ymax=max(0,fr[:,1].max())
lim=max(xmax-xmin,ymax-ymin)*1.15
ax.set(xlabel='Along observed motion (km)',ylabel='Perpendicular motion (km)')
ax.set_aspect('equal',adjustable='box'); ax.set_xticks([0,200,400]); ax.set_yticks([-100,0,100,200])
ax.set_xlim(xmin-.05*lim,xmin+.95*lim); ax.set_ylim(ymin-.05*lim,ymin+.95*lim)
ax.legend(handles=[Line2D([],[],color=C['full'],marker='s',mfc='white',ms=4.5,label='Full'),handle('early','Early half'),Line2D([],[],color='#222222',marker='+',label='Observed')],
    loc='upper right',bbox_to_anchor=(1,.85),fontsize=7.5,handlelength=1.2,labelspacing=.25,borderpad=.15)
ax.plot([fm[0],fm[0]],[fm[1],90],color='#999999',ls=':',lw=.55)
ax.plot([lo,lo],[0,90],color='#999999',ls=':',lw=.55)
ax.annotate('',xy=(lo,90),xytext=(fm[0],90),arrowprops={'arrowstyle':'<->','lw':.7,'color':'#444444'})
ax.text((lo+fm[0])/2,105,f'{lo-fm[0]:.0f} km',ha='center',va='bottom',fontsize=7.5)
# Preserve the true endpoints. A larger hollow square exposes the full mean
# behind the early-half diamond; the inset resolves their actual small separation.
ax.plot(*fm,marker='s',ms=6.7,mfc='none',mec=C['full'],mew=1.0,ls='none',zorder=8)
detail=ax.inset_axes([.12,.55,.30,.30])
detail.set_facecolor('white')
detail.axhline(0,color='#DDDDDD',lw=.5,zorder=0)
detail.axvline(0,color='#DDDDDD',lw=.5,zorder=0)
delta=sm-fm
detail.plot(0,0,marker='s',ms=4.8,mfc='white',mec=C['full'],mew=.9,ls='none')
detail.plot(*delta,marker='D',ms=3.4,color=C['early'],ls='none')
detail.set(xlim=(-.8,1.2),ylim=(-1.5,.5),xticks=[0,1],yticks=[-1,0])
detail.set_aspect('equal',adjustable='box')
detail.tick_params(length=2,pad=1,labelsize=7.5)
detail.text(.97,.97,'km',transform=detail.transAxes,ha='right',va='top',fontsize=7.5)
detail.annotate('',xy=delta,xytext=(0,0),arrowprops={'arrowstyle':'->','color':'#888888','lw':.6,'shrinkA':4,'shrinkB':4})
ax.annotate('',xy=fm,xycoords='data',xytext=(.43,.55),textcoords='axes fraction',arrowprops={'arrowstyle':'-','color':'#888888','lw':.55,'linestyle':':'},zorder=2)
(O/'endpoint_detail.json').write_text(json.dumps({'reference':'full-ensemble endpoint','full_endpoint_km':fm.tolist(),'early_endpoint_km':sm.tolist(),'early_minus_full_km':delta.tolist(),'display_displacement_added':False,'same_Hagibis_case':True},indent=2))
(O/'segment_geometry.json').write_text(json.dumps({'case':'Hagibis 2019, 2019-10-07 12UTC','segment':[108,120],
    'member_along_min':float(fr[:,0].min()),'member_along_max':float(fr[:,0].max()),'members':len(fr),'mean_member_length':float(L),'mean_vector_length':float(s),'observed_length':float(lo),
    'mean_vector_observed_frame':fm.tolist(),'early_vector_observed_frame':sm.tolist(),'identity_residual':(terms.sum(0)-(fm-[lo,0])).tolist()},indent=2))

termnames=['member-length term','cancellation term','direction term','total ATE']
for i,center in enumerate(['ECMWF','NCEP']):
    ax=axs[1+i];title(ax,chr(98+i),center+' · lag components')
    ax.axvline(0,color='#888888',lw=.7);ax.axhline(2.5,color='#DDDDDD',lw=.6)
    for j,term in enumerate(termnames):
        for kind,off,key in [('Full mean',-.13,'full'),('Selection change',.13,'early')]:
            q=row(parts,center=center,kind=kind,term=term)
            ax.errorbar(q['mean'],j+off,xerr=[[q['mean']-q.ci_low],[q.ci_high-q['mean']]],fmt=M[key],color=C[key],ms=3.5,capsize=2,lw=1,zorder=3)
        if j==3:
            q=row(parts,center=center,kind='Selection change',term=term)
            ax.annotate(f'{q["mean"]:+.1f}',(q['mean'],j+.13),xytext=(-25,-10),textcoords='offset points',color=C['early'],ha='center',fontsize=8)
    ax.set(yticks=range(4),yticklabels=['Length','Averaging','Direction','Total'],ylim=(4.9,-.6),xlim=(-225,32),xticks=[-200,-100,0],xlabel='Contribution or change (km)')
    ax.tick_params(axis='y',length=0,pad=2)
    ax.legend(handles=[handle('full','Full mean'),handle('early','Selected − full')],loc='lower left',fontsize=7.5,ncol=1,handlelength=1.2,labelspacing=.25,borderpad=.1)


fig.text(.5,.025,'L: mean member length; s: length of the mean vector; obs: observed displacement length',ha='center',fontsize=8)
save(fig,axs,'evidence_displacement_budget')

