from pathlib import Path
import json,hashlib
from plot_fonts import configure_font
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt,font_manager
from matplotlib.patches import Circle,Ellipse,Rectangle,FancyArrowPatch
from matplotlib.lines import Line2D
from audit_panel_alignment import require_matplotlib_panel_alignment
import os
B=Path(__file__).resolve().parent.parent;O=Path(os.environ['REPRO_OUTPUT'])/'fig1';O.mkdir(parents=True,exist_ok=False)
FONT=configure_font()
plt.rcParams.update({'font.family':FONT,'font.size':8,'pdf.fonttype':42,'svg.fonttype':'none','mathtext.fontset':'stix','figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white','savefig.transparent':False,'legend.frameon':False,'axes.linewidth':.65})
C='#178F9F';S='#D48732';P='#8062A6';G='#B9BEC2';K='#222222'
width_mm=183
W,H=180,108
fig=plt.figure(figsize=(W/25.4,H/25.4))
def axmm(x,y,w,h):return fig.add_axes([x/W,y/H,w/W,h/H])
a=axmm(8,54,58,44);b=axmm(77,54,95,44);c=axmm(8,4,164,36)
axes=[a,b,c]
def title(ax,l,t):
 r=ax.get_position();ax.text(0,1+2/(r.height*H),l,transform=ax.transAxes,fontsize=10,weight='bold',va='bottom')
 ax.text(4.5/(r.width*W),1+2/(r.height*H),t,transform=ax.transAxes,fontsize=9.5,va='bottom')
def curve(start,control,end):
 t=np.linspace(0,1,100)[:,None];return (1-t)**2*start+2*t*(1-t)*control+t*t*end
def cyclone(ax,xy,r=.10):
 ax.add_patch(Circle(xy,r,facecolor='white',edgecolor=K,lw=1,zorder=8))
 ang=np.pi/4;rot=np.array([[np.cos(ang),-np.sin(ang)],[np.sin(ang),np.cos(ang)]])
 for arm in [np.array([[0,r],[.24,r]]),np.array([[-.24,-r],[0,-r]])]:
  v=arm@rot.T+xy;ax.plot(*v.T,color=K,lw=1,zorder=8)
rawerr=np.array([[-.22,.18],[.16,.25],[-.10,-.24],[.27,-.10],[-.58,.62],[.68,.35],[.60,-.48],[-.46,-.59]])
keep=np.argsort(np.linalg.norm(rawerr,axis=1),kind='stable')[:4]
assert set(keep)=={0,1,2,3}
anchor=np.array([2.15,1.1]);init=np.array([.25,.05]);raw=anchor+rawerr
colors=[S if i in keep else G for i in range(8)]
title(a,'a','Select using raw 24-h errors')
a.set(xlim=(-.05,3.15),ylim=(-.35,2.63));a.set_aspect('equal',adjustable='datalim');a.axis('off')
for i in range(8):
 past=curve(init,(init+raw[i])/2+np.array([0,-.18]),raw[i])
 a.plot(*past.T,color=colors[i],lw=.75 if i in keep else .55,alpha=.8 if i in keep else .6,zorder=1)
 a.plot(*raw[i],marker='o',ms=4.4,mfc='white',mec=colors[i],mew=1,zorder=4)
obs=curve(init,np.array([1.2,.18]),anchor);a.plot(*obs.T,color=K,lw=1.2,zorder=3)
cyclone(a,init)
a.plot(*anchor,marker='s',ms=4.8,mfc='white',mec=K,mew=1,zorder=6)
a.text(.67,-.18,'0 h',fontsize=7.8,ha='center',va='center')
a.annotate('Observed\n24-h position',xy=anchor,xytext=(1.35,2.1),fontsize=8,va='top',linespacing=1.2,arrowprops={'arrowstyle':'-','color':K,'lw':.6,'connectionstyle':'arc3,rad=0'})
a.annotate('Retained half',xy=raw[0],xytext=(.05,2.05),fontsize=8,color=S,ha='left',va='top',arrowprops={'arrowstyle':'-','color':S,'lw':.65})
fig.text(8/W,47/H,'Hollow circles: raw member positions',fontsize=7.5,ha='left')

# Later paths are the member displacements after each own 24-h error is removed.
title(b,'b','Translate each member; evaluate later')
b.set(xlim=(-.45,8.63),ylim=(-.6,4.4));b.set_aspect('equal',adjustable='datalim');b.axis('off')
ends=np.array([[3.0,2.7],[3.4,3.1],[2.7,3.2],[3.1,2.6],[3.8,3.9],[1.6,1.1],[3.9,1.8],[1.5,2.4]])
start=np.zeros(2);truth=np.array([3.7,3.8]);controls=np.array([1.35,.45])+.20*(ends-ends.mean(0))
paths=np.array([curve(start,controls[i],ends[i]) for i in range(8)])
for i in range(8):
 b.plot(*paths[i].T,color=colors[i],lw=.7 if i in keep else .55,alpha=.68 if i in keep else .65,zorder=1)
 b.plot(*ends[i],marker='o',ms=3.7,mfc=colors[i],mec=colors[i],mew=.3,alpha=.85,zorder=3)
obs=curve(start,np.array([1.45,.50]),truth);b.plot(*obs.T,color=K,lw=1.3,zorder=2)
fullpath=paths.mean(0);selpath=paths[keep].mean(0)
b.plot(*fullpath.T,color=C,lw=1.5,zorder=4)
b.plot(*selpath.T,color=S,lw=1.5,ls=(0,(4,2)),zorder=5)
full=ends.mean(0);sel=ends[keep].mean(0)
b.plot(*full,marker='s',ms=5,mfc=C,mec='white',mew=.5,zorder=7)
b.plot(*sel,marker='D',ms=5,mfc=S,mec='white',mew=.5,zorder=8)
b.plot(*start,marker='s',ms=4.8,mfc='white',mec=K,mew=1,zorder=8)
b.plot(*truth,marker='s',ms=4.8,mfc='white',mec=K,mew=1,zorder=8)
b.text(-.18,-.35,'Observed 24-h anchor',fontsize=7.8,ha='left',va='top')
for text,xy,pos,col in [('Observed 120 h',truth,(5.0,4.0),K),('Selected mean',sel,(5.0,3.1),S),('Full mean',full,(5.0,2.2),C)]:
 b.annotate(text,xy=xy,xytext=pos,fontsize=8,color=col,ha='left',va='center',arrowprops={'arrowstyle':'-','color':col,'lw':.65,'shrinkA':3,'shrinkB':4})
mid=(sel+truth)/2;off=np.array([-.48,.30])
for endpoint in [sel,truth]:b.plot([endpoint[0],endpoint[0]+off[0]],[endpoint[1],endpoint[1]+off[1]],color='#777777',ls=':',lw=.5,zorder=2)
b.annotate('',xy=truth+off,xytext=sel+off,arrowprops={'arrowstyle':'<->','color':K,'lw':.7,'shrinkA':1,'shrinkB':1},zorder=7)
b.text(.55,3.65,'Later position\nerror remains',fontsize=7.8,va='center',linespacing=1.15)
b.text(5.0,.85,'Orange: retained members\nGray: other members',fontsize=7.8,linespacing=1.4,va='center')
fig.text(77/W,48/H,"Subtract each member's own 24-h position error.",fontsize=8.3,ha='left')
# Visual step arrow between spatial panels, separate from geographic displacement.
fig.add_artist(FancyArrowPatch((68/W,76/H),(73/W,76/H),transform=fig.transFigure,arrowstyle='->',mutation_scale=8,lw=.9,color=K))

title(c,'c','Three complementary comparisons')
c.set(xlim=(0,167),ylim=(0,40));c.axis('off')
for x,col,head in [(0,C,'Realized mean-track change'),(57,S,'Order and weight shape'),(114,P,'Member-range coverage')]:
 c.text(x,37,head,fontsize=9.5,color=col,va='top')
for label,xy,col,marker in [('Full',(12,29),C,'s'),('Half',(12,12),S,'D')]:
 c.text(0,xy[1],label,fontsize=8,va='center');c.plot([xy[0],42],[xy[1],21],color=col,lw=1.1);c.plot(*xy,marker=marker,ms=5,color=col)
c.plot(42,21,marker='+',ms=7,color=K);c.text(44,21,'Obs.',fontsize=8,va='center')
for j,(label,order) in enumerate([('Early',[1,2,3,4]),('Verifying',[3,1,4,2]),('Random',[4,2,1,3])]):
 y=28-j*4
 c.text(57,y,label,fontsize=8.5,va='center')
 for k,member in enumerate(order):c.text(79+k*6,y,str(member),fontsize=8.5,ha='center',va='center')
weights=[np.array([.4,.4,.1,.1]),np.array([.4,.3,.2,.1])]
for x,w,label in zip([58,84],weights,['Two-level','Monotone']):
 assert np.all(w>=0) and np.isclose(w.sum(),1)
 c.step(x+np.arange(5)*4,9+np.r_[w,w[-1]]*12,where='post',color=S,lw=1)
 c.text(x,7,label,fontsize=8.5,va='center')
c.text(57,3,'Either direction; schematic',fontsize=8.5,va='center')
center=(137,20);c.add_patch(Circle(center,12,fill=False,ec=P,lw=1));c.add_patch(Ellipse(center,35,19,angle=20,fill=False,ec=P,lw=1,ls='--'))
c.plot(*center,marker='s',ms=4,mfc='white',mec=K);c.plot(159,29,marker='+',ms=7,color=K);c.text(161,29,'Obs.',fontsize=8,va='center')
for x,s in [(0,'24-h update'),(57,'Future-truth optimization'),(114,'Equal-size random reference')]:c.text(x,-2,s,fontsize=8.5,color='#555555',va='center')
fig.canvas.draw()
# Semantic checks for the schematic construction, not empirical claims.
assert np.allclose(fullpath[-1],full) and np.allclose(selpath[-1],sel)
assert all(np.allclose(p[0],start) for p in paths)
assert min(t.get_fontsize() for ax in axes for t in ax.texts)>=7.5
require_matplotlib_panel_alignment(fig,json_out=O/'schematic_alignment.json',axes=axes,panel_ids=list('abc'),row_groups=[['a','b']],column_groups=[],exemptions=[{'panels':['a','b'],'checks':['panel-width'],'reason':'The later-track panel is intentionally wider for directly attached labels.'}],tolerance_pt=1.5,gutter_tolerance_pt=1.5,require_panel_labels=True,strict=True)
fig.savefig(O/'Fig01_StudyDesign.pdf',dpi=600)
fig.savefig(O/'Fig01_StudyDesign.svg',dpi=600)
fig.savefig(O/'Fig01_StudyDesign.png',dpi=600)
(O/'schematic_semantics.json').write_text(json.dumps({'status':'PASS','construction':'Original idealized vector trajectories, not observed or typical. No physical mechanism inferred.','raw_errors':rawerr.tolist(),'selected':keep.tolist(),'later_endpoints':ends.tolist(),'full_mean':full.tolist(),'selected_mean':sel.tolist(),'memberwise_translation':'Each original raw 24-h offset is subtracted from that entire member path; diagram b displays the resulting displacements.','display_order':'Raw 24-h ranking before translation. Membership stays fixed.','weighting':'Later truth is used by every optimized weighting family; verifying order is not forecast-time information.','coverage':'MW/CT mean identity restricted to each fixed equal-weight group. Range definitions are not calibrated probabilities.','panel_c':'Illustrative error vectors, rank permutations and equal-center region boundaries. Three complementary tests, not a sequential workflow. No actual effect size is encoded.', 'new_empirical_analysis':False,'matplotlib':matplotlib.__version__},indent=2))
print('Schematic rendered, semantic checks and alignment PASS',flush=True)
