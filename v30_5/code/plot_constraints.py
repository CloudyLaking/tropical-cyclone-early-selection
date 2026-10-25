"""Derived-data reconstruction: paired real examples, storm distributions and random-order controls."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import argparse,json,hashlib
from plot_fonts import configure_font
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Polygon,Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from audit_panel_alignment import require_matplotlib_panel_alignment

HERE=Path(__file__).resolve().parent
SHA=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
COL=dict(early='#E98B2A',verifying='#8B5FBF',free='#222222',full='#0095A6',random='#777777')
FONT=configure_font()
plt.rcParams.update({'font.family':FONT,'font.size':7.5,'axes.labelsize':7.5,'axes.titlesize':8,
 'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':7.5,'axes.linewidth':.65,
 'xtick.major.width':.65,'ytick.major.width':.65,'xtick.major.size':2.5,'ytick.major.size':2.5,
 'axes.spines.top':True,'axes.spines.right':True,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
 'savefig.facecolor':'white','figure.facecolor':'white'})

def title(ax,letter,label):
    ax.text(0,1.045,letter,transform=ax.transAxes,fontsize=8,fontweight='bold',va='bottom')
    ax.text(.075,1.045,label,transform=ax.transAxes,fontsize=8,va='bottom')

def main(data,out):
    out.mkdir(parents=True,exist_ok=True)
    assert json.loads((data/'data_receipt.json').read_text())['status']=='PASS'
    read=lambda n:pd.read_csv(data/n,float_precision='round_trip')
    members,sets,points,storms,summary,random=[read(n) for n in ['example_members.csv','example_sets.csv','example_points.csv','storm_contrasts.csv','contrast_summary.csv','random_contrasts.csv']]
    width_mm=183
    fig=plt.figure(figsize=(width_mm/25.4,140/25.4))
    axes=[fig.add_axes([x/183,y/140,w/183,39/140]) for x,w in [(13,42),(82,52),(144,35)] for y in [85,28]]
    fig.text(.5,137/140,'a,d: Hagibis · 7 Oct 2019, 12 UTC     b,c,e,f: all complete-cohort forecasts',ha='center',va='center',fontsize=8)
    fig.text(.5,132/140,'All optima use later truth; M permits increasing or decreasing weights',ha='center',va='center',fontsize=7.5)
    contrasts=['half_minus_share','share_minus_split','split_minus_levels','early_minus_verifying','verifying_minus_free']
    labels=['Early half −\nfixed split','Fixed split −\noptimal split','Optimal split −\nearly M','Early M −\nverifying M','Verifying M −\nfree']
    colors=[COL['early']]*3+[COL['verifying'],COL['free']]
    limits=[storms[contrasts].min().min(),storms[contrasts].max().max()]
    span=limits[1]-limits[0];gain_lim=(min(-10,limits[0]-.06*span),limits[1]+.08*span)
    ymin=float(random.ci_low.min());ymax=float(random.ci_high.max());dy=ymax-ymin
    for c,center in enumerate(['ECMWF','NCEP']):
        top,mid,bot=axes[c],axes[2+c],axes[4+c]
        m=members[members.center==center];s=sets[sets.center==center];p=points[points.center==center].set_index('family')
        for family,color,style,alpha in [('H',COL['free'],'-',.025),('early_D',COL['early'],'--',.10),('verifying_D',COL['verifying'],'-',.10)]:
            q=s[(s.family==family)&(s.hull_order>=0)].sort_values('hull_order')[['x','y']].to_numpy()
            if len(q)>2:
                top.add_patch(Polygon(q,closed=True,facecolor=color,alpha=alpha,edgecolor='none',zorder=1));q=np.vstack([q,q[0]])
            top.plot(q[:,0],q[:,1],color=color,lw=.8,ls=style,zorder=2)
        top.scatter(m.x,m.y,s=8,c='#9BA3A7',alpha=.70,edgecolors='white',linewidths=.2,zorder=3)
        for family,color,marker,size,face in [('Full',COL['full'],'s',4.2,COL['full']),('early_D',COL['early'],'s',4.6,COL['early']),('verifying_D',COL['verifying'],'o',5.7,'white'),('H',COL['free'],'D',3.0,COL['free'])]:
            q=p.loc[family];top.plot(q.x,q.y,ls='none',marker=marker,ms=size,color=color,mfc=face,mew=.85,zorder=5)
        top.plot(0,0,marker='+',ms=6,mew=1,color=COL['free'],zorder=6)
        xy=np.vstack([members[['x','y']],np.zeros((1,2))]);extent=xy.max(0)-xy.min(0);centerxy=(xy.max(0)+xy.min(0))/2
        ratio=42/39;hy=max(extent[1]/2,extent[0]/2/ratio)*1.15
        top.set_xlim(centerxy[0]-ratio*hy,centerxy[0]+ratio*hy);top.set_ylim(centerxy[1]-hy,centerxy[1]+hy)
        top.set_aspect('equal',adjustable='box');top.xaxis.set_major_locator(MaxNLocator(3));top.yaxis.set_major_locator(MaxNLocator(4))
        top.axhline(0,color='#dddddd',lw=.4,zorder=0);top.axvline(0,color='#dddddd',lw=.4,zorder=0)
        top.set_xlabel('Projected eastward error (km)',labelpad=3);top.set_ylabel('Projected northward error (km)',labelpad=3)
        title(top,'ad'[c],f'{center} · decreasing D')
        st=storms[storms.center==center].sort_values('sid');su=summary[summary.center==center].set_index('contrast')
        for i,(key,color) in enumerate(zip(contrasts,colors)):
            y=[4,3,2,1,0][i];q=su.loc[key];values=st[key].to_numpy()
            b=dict(med=q.q50,q1=q.q25,q3=q.q75,whislo=q.q05,whishi=q.q95,fliers=[])
            mid.bxp([b],positions=[y],widths=.36,vert=False,showfliers=False,patch_artist=True,
                boxprops=dict(facecolor=color,edgecolor=color,alpha=.20,linewidth=.8),medianprops=dict(color=color,lw=1),
                whiskerprops=dict(color=color,lw=.8),capprops=dict(color=color,lw=.8))
            jitter=np.array([(int(hashlib.sha256((sid+key).encode()).hexdigest()[:8],16)/4294967295-.5)*.27 for sid in st.sid])
            mid.scatter(values,y+jitter,s=5,color=color,alpha=.38,edgecolors='none',zorder=3)
            mid.errorbar(q['mean'],y,xerr=[[q['mean']-q.ci_low],[q.ci_high-q['mean']]],fmt='D',ms=3.4,
                color=COL['free'],mfc='white',mew=.75,capsize=1.4,lw=.9,zorder=5)
            mid.text(.98,y+.38,f"{q['mean']:.1f} [{q.ci_low:.1f}, {q.ci_high:.1f}]",transform=mid.get_yaxis_transform(),ha='right',va='baseline',fontsize=7.5,color=COL['free'])
        mid.axvline(0,color='#9a9a9a',ls='--',lw=.65,zorder=0)
        mid.set_yticks([0,1,2,3,4],labels[::-1]);mid.set_ylim(-.55,4.85);mid.set_xlim(*gain_lim);mid.xaxis.set_major_locator(MaxNLocator(4))
        mid.set_xlabel('Paired distance reduction (km)',labelpad=3)
        title(mid,'be'[c],'Weight / order contrasts')
        for metric,color,mark,offset,ls in [('D_early_minus_random',COL['early'],'s',-.9,'-'),('M_early_minus_random',COL['free'],'o',.9,'--')]:
            q=random[(random.center==center)&(random.metric==metric)].sort_values('lead')
            bot.errorbar(q.lead+offset,q['mean'],yerr=[q['mean']-q.ci_low,q.ci_high-q['mean']],fmt=mark,
                ls=ls,color=color,mfc='white' if mark=='o' else color,ms=3.6,mew=.75,capsize=1.8,lw=.9)
        bot.axhline(0,color='#9a9a9a',ls='--',lw=.7,zorder=0)
        bot.set_xlim(43,125);bot.set_xticks([48,72,96,120]);bot.set_ylim(ymin-.12*dy,ymax+.14*dy);bot.yaxis.set_major_locator(MaxNLocator(4))
        bot.set_xlabel('Lead time (h)',labelpad=3)
        title(bot,'cf'[c],'Early − random (km)')
        for ax in [top,mid,bot]:
            ax.tick_params(direction='out',pad=2)
            for sp in ax.spines.values():sp.set_visible(True)
    geom=[Line2D([],[],color='#9BA3A7',marker='o',ls='none',ms=3,label='Members'),Line2D([],[],color=COL['full'],marker='s',ls='none',ms=4,label='Full equal mean'),Line2D([],[],color=COL['free'],marker='+',ls='none',ms=5,label='Truth'),Line2D([],[],color=COL['early'],ls='--',marker='s',ms=3.5,label='Early-order decreasing (D)'),Line2D([],[],color=COL['verifying'],marker='o',mfc='white',ms=4,label='Verifying-order decreasing (D)'),Line2D([],[],color=COL['free'],marker='D',ms=3,label='Free weights')]
    fig.legend(handles=geom,loc='center',bbox_to_anchor=(.5,12/140),ncol=3,frameon=False,columnspacing=1.3,handletextpad=.5,handlelength=1.6,labelspacing=.4)
    dist=[Line2D([],[],color='#777777',marker='o',ls='none',ms=2.8,label='All storms'),Patch(fc='#dddddd',ec='#777777',label='Quartiles; 5–95% whiskers'),Line2D([],[],color=COL['free'],marker='D',mfc='white',ms=3.5,label='Mean; paired 95% CI')]
    # Box, point and interval conventions are stated next to their values in the caption.
    control=[Line2D([],[],color=COL['early'],marker='s',ms=3.5,label='Decreasing only (D)'),Line2D([],[],color=COL['free'],ls='--',marker='o',mfc='white',ms=3.5,label='Either direction (M)')]
    fig.legend(handles=control,loc='lower center',bbox_to_anchor=(.5,0/140),ncol=2,frameon=False,columnspacing=2.,handletextpad=.5,handlelength=1.8)
    fig.canvas.draw()
    audit=require_matplotlib_panel_alignment(fig,axes=axes,panel_ids=list('adbecf'),row_groups=[['a','b','c'],['d','e','f']],column_groups=[['a','d'],['b','e'],['c','f']],exemptions=[dict(panels=list('abcdef'),checks=['panel-width','horizontal-gutter'],reason='Explicit unequal columns: 42 mm geometric examples, 52 mm distribution panels, 35 mm leadwise contrasts; gaps reserve room for paired-rule labels. Within-column widths and all row edges retain the 1.5pt tolerance.')],json_out=out/'Fig03_Constraints.alignment.json',overlay_svg=out/'Fig03_Constraints.alignment.svg',tolerance_pt=1.5,require_panel_labels=True)
    fig.savefig(out/'Fig03_Constraints.pdf',dpi=600)
    fig.savefig(out/'Fig03_Constraints.svg',dpi=600)
    fig.savefig(out/'Fig03_Constraints.png',dpi=600)
    caption='Weight rules and ordering information. Rows show ECMWF/NCEP. (a,d) Actual 120-h member errors for Hagibis (2019), initialized 12 UTC 7 October 2019. Boundaries contain means attainable with free nonnegative weights or decreasing weights (D) under early/verifying order; markers show their closest points to truth. These are feasible means, not probability regions. (b,e) Paired distance reductions across 64/100 storms. Orange contrasts hold early order fixed while optimizing shares, dividing rank and monotone levels; purple replaces early order with verifying order, and black removes ordering constraints. Fixed split and Optimal split permit either group to receive greater per-member weight; M permits either monotone direction, whereas (a,d) show D only. Positive values favor the second named rule. Boxes/points show storm distributions; diamonds and bars show means and paired 95\\% intervals, printed at right. The contrasts are not additive causal components. (c,f) Early-minus-random optimal distances across leads, with original paired 95\\% intervals; negative favors early order. Every optimized weight uses later truth independently at each target. Random results average 128 separately optimized permutations, never the best permutation.'
    (out/'Fig03_Constraints.caption.md').write_text(caption,encoding='utf-8')
    receipt=dict(status='RENDERED',alignment=audit['verdict'],width_mm=183,height_mm=140,font='Arial',minimum_font_pt=7.5,script_sha256=SHA(__file__),data_receipt_sha256=SHA(data/'data_receipt.json'),new_optimization=False,all_storms_plotted=True,case_selection='previously paired Hagibis, 2019-10-07 12UTC')
    (out/'render_receipt.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2));plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();main(a.data,a.out)
