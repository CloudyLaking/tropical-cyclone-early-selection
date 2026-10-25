"""SERVER ONLY: two independent manuscript figures from frozen ranking summaries.

Invoke separately with --figure 1 or 2. Source data remain read-only.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from matplotlib.ticker import MaxNLocator,FixedLocator
from audit_panel_alignment import require_matplotlib_panel_alignment

HERE=Path(__file__).resolve().parent
DATA=HERE.parent/'data/ranking'
COL=dict(early='#E98B2A',late='#8B5FBF',random='#777777',free='#222222',full='#0095A6')
SHA=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
width_mm=183
from plot_fonts import configure_font
FONT=configure_font()
plt.rcParams.update({'font.family':FONT,'font.size':7.5,'axes.labelsize':7.5,
    'axes.titlesize':8,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':7.5,
    'axes.linewidth':.65,'xtick.major.width':.65,'ytick.major.width':.65,
    'xtick.major.size':2.5,'ytick.major.size':2.5,'axes.unicode_minus':True,
    'axes.spines.top':True,'axes.spines.right':True,'pdf.fonttype':42,'ps.fonttype':42,
    'svg.fonttype':'none','savefig.facecolor':'white','figure.facecolor':'white'})

def title(ax,letter,label):
    ax.text(0,1.045,letter,transform=ax.transAxes,fontsize=8,fontweight='bold',va='bottom')
    ax.text(.065,1.045,label,transform=ax.transAxes,fontsize=8,va='bottom')

def neat(ax):
    ax.tick_params(direction='out',pad=2)
    for spine in ax.spines.values():spine.set_visible(True)

def export(fig,out,name,axes,rows,columns,source_files,caption):
    out.mkdir(parents=True,exist_ok=True)
    fig.canvas.draw()
    report=require_matplotlib_panel_alignment(fig,axes=axes,panel_ids=list('abcdef')[:len(axes)],
        row_groups=rows,column_groups=columns,json_out=out/(name+'.alignment.json'),
        overlay_svg=out/(name+'.alignment.svg'),tolerance_pt=1.5,require_panel_labels=True)
    fig.savefig(out/(name+'.pdf'),dpi=600)
    fig.savefig(out/(name+'.svg'),dpi=600)
    fig.savefig(out/(name+'.png'),dpi=600)
    (out/(name+'.caption.md')).write_text(caption,encoding='utf-8')
    result=dict(status='RENDERED',figure=name,
        scientific_scope='Frozen summaries; terminology revision only',width_mm=fig.get_figwidth()*25.4,
        height_mm=fig.get_figheight()*25.4,font=FONT,font_min_pt=7.5,
        source_sha256={str(p.relative_to(HERE.parent)):SHA(p) for p in source_files},script_sha256=SHA(__file__),
        alignment_verdict=report['verdict'],no_new_statistics=True,
        review_notes=['All case members retained. No samples selected from V28 outcome.',
                      'All confidence intervals are precomputed pointwise exploratory 95% paired storm intervals.'])
    (out/(name+'.provenance.json')).write_text(json.dumps(result,indent=2))
    plt.close(fig)

def getrow(df,center,lead,metric):
    q=df[(df.center==center)&(df.lead.astype(str)==str(lead))&(df.metric==metric)]
    assert len(q)==1,(center,lead,metric)
    return q.iloc[0]

def figure2(out):
    files=[DATA/'summary.csv',HERE.parent/'data/si/legacy_group.csv',HERE.parent/'data/si/legacy_monotone.csv']
    new,group,old=[pd.read_csv(p,float_precision='round_trip') for p in files]
    clean=out/'source_data';clean.mkdir(parents=True,exist_ok=True);used=[]
    fig=plt.figure(figsize=(183/25.4,166/25.4))
    axes=[fig.add_axes([x,y,.325,h]) for y,h in [(.705,.255),(.3925,.215),(.145,.15)] for x in [.235,.66]]
    method_rows=[('Full mean',group,'full',COL['full'],'s'),
        ('Early half',group,'actual',COL['early'],'^'),
        ('Fixed split',group,'group',COL['early'],'v'),
        ('Best split',old,'two_level',COL['early'],'D'),
        ('Early M',new,'M_early',COL['early'],'s'),
        ('Verifying M',new,'M_late',COL['late'],'o'),
        ('Free nonnegative',new,'H',COL['free'],'D')]
    contrasts=[('Best split − early M',old,'multi_level_gain',COL['early']),
        ('Early M − verifying M',new,'M_early_minus_late',COL['late']),
        ('Verifying M − free',new,'M_late_minus_H',COL['free']),
        ('Early M − random M',new,'M_early_minus_random',COL['random'])]
    for c,center in enumerate(['ECMWF','NCEP']):
        top,mid,bot=axes[c],axes[2+c],axes[4+c]
        for i,(label,df,key,color,mark) in enumerate(method_rows):
            row=getrow(df,center,120,key);y=6-i
            top.errorbar(row['mean'],y,xerr=[[row['mean']-row.ci_low],[row.ci_high-row['mean']]],fmt=mark,
                         color=color,ms=4,capsize=2,lw=.9,mew=.8)
            used.append(dict(panel='ab'[c],display=label,source_metric=key,**row.to_dict()))
        top.set_yticks(range(7),[r[0] for r in method_rows][::-1] if c==0 else ['']*7)
        top.set_ylim(-.6,6.6);top.set_xlim(0,500);top.set_xticks([0,250,500])
        top.set_xlabel('Mean projected distance (km)',labelpad=3);title(top,'ab'[c],center+' · 120 h')
        for label,key,color,marker,offset,style in [('Early M','M_early',COL['early'],'s',-1.0,'-'),
              ('Verifying M','M_late',COL['late'],'o',-.33,'-'),
              ('Random M','M_random',COL['random'],'D',.33,'--'),
              ('Free nonnegative','H',COL['free'],'^',1.0,':')]:
            q=pd.DataFrame([getrow(new,center,h,key) for h in [48,72,96,120]])
            mid.errorbar(np.array([48,72,96,120])+offset,q['mean'],yerr=[q['mean']-q.ci_low,q.ci_high-q['mean']],
                 fmt=marker,ls=style,ms=3.5,lw=.9,color=color,capsize=1.4,mew=.65,mfc='white' if key=='M_random' else color)
            for _,row in q.iterrows():used.append(dict(panel='cd'[c],display=label,source_metric=key,**row.to_dict()))
        mid.set_xlim(43,125);mid.set_xticks([48,72,96,120]);mid.set_ylim(0,280);mid.set_yticks([0,100,200])
        mid.set_xlabel('Lead time (h)',labelpad=3);mid.set_ylabel('Mean distance (km)',labelpad=3)
        title(mid,'cd'[c],center+' · change with lead time')
        for i,(label,df,key,color) in enumerate(contrasts):
            row=getrow(df,center,120,key);y=3-i
            bot.errorbar(row['mean'],y,xerr=[[row['mean']-row.ci_low],[row.ci_high-row['mean']]],fmt='o',
                         color=color,ms=4,capsize=2,lw=.9,mew=.8)
            used.append(dict(panel='ef'[c],display=label,source_metric=key,**row.to_dict()))
        bot.axvline(0,color='#aaaaaa',ls='--',lw=.7,zorder=0)
        bot.set_yticks(range(4),[r[0] for r in contrasts][::-1] if c==0 else ['']*4)
        bot.set_ylim(-.6,3.6);bot.set_xlim(-35,170);bot.set_xticks([0,80,160])
        bot.set_xlabel('Paired distance difference (km)',labelpad=3);title(bot,'ef'[c],center+' · paired differences at 120 h')
        for ax in [top,mid,bot]:neat(ax)
    handles=[Line2D([],[],color=COL['early'],marker='s',ms=4,label='Early M'),
             Line2D([],[],color=COL['late'],marker='o',ms=4,label='Verifying M'),
             Line2D([],[],color=COL['random'],ls='--',marker='D',mfc='white',ms=4,label='Random M'),
             Line2D([],[],color=COL['free'],ls=':',marker='^',ms=4,label='Free nonnegative')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.56,.012),ncol=2,frameon=False,
               columnspacing=1.6,handlelength=2.0,handletextpad=.6,labelspacing=.45)
    pd.DataFrame(used).to_csv(clean/'plotted_summary_rows.csv',index=False)
    caption='''Figure 2. Attainable track distances under early, verifying and random ranking constraints.

(a,b) Storm-equal 120-h mean projected distances for the same 806 initializations (ECMWF: 300 initializations/64 storms; NCEP: 506/100). Full and early-half equal means are fixed products. The fixed-half comparison optimizes the weight shared between the two equal-weight half means; the two-level family additionally chooses the split along the early order. Early- and verifying-order monotone families (M) optimize arbitrary nonnegative normalized weights constrained to either nonincreasing or nondecreasing order, selecting the better direction for each initialization and lead. Free nonnegative use the full member convex hull. All optimizations use future verifying truth, independently at each lead; they are diagnostic attainability comparisons, not operational forecasts.

(c,d) M-family and unconstrained distances across 48–120 h. Random-order M averages the optimized result of each of 128 independently generated full orders within initialization; it does not choose the best order. The same 128 orders are retained across all four leads. Small horizontal offsets separate error bars. (e,f) Saved paired contrasts at 120 h; the first is from the preceding two-level/monotone comparison, and the remaining three use the new controlled ranking comparison. Positive distance differences favor the second named method. Early versus verifying sets are not nested, so their difference is not a causal attribution. Points are storm-equal means; all bars are pointwise exploratory 95% percentile intervals from the same 5000 paired whole-storm draws. No multiplicity correction. D-family results and all non-significant comparisons remain in the tables. Scientific independent acceptance of V02 is pending in this initial rendering.
'''
    export(fig,out,'Fig02_RankingPopulation',axes,[['a','b'],['c','d'],['e','f']],
           [['a','c','e'],['b','d','f']],files,caption)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--figure',choices=['2'],required=True);ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    figure2(args.out)
    print(json.dumps(dict(status='rendered',figure=args.figure,out=str(args.out))))
