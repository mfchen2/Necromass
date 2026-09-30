from pathlib import Path
import importlib.util
import pandas as pd
import numpy as np
import textwrap
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

OUT = Path(__file__).resolve().parent
SOURCE = OUT.parent / 'manual-20260609-a8/presentations/mofa_composite_abcdef/assets/reproduce_panel.py'
spec = importlib.util.spec_from_file_location('original', SOURCE)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
plt.rcParams.update({'font.size':14, 'axes.labelsize':16, 'xtick.labelsize':14,
                     'ytick.labelsize':14, 'pdf.fonttype':42, 'svg.fonttype':'none'})
fig, axes = plt.subplots(3,2,figsize=(26,29),gridspec_kw={'height_ratios':[1,1,2]})
fig.subplots_adjust(left=.13,right=.95,top=.95,bottom=.08,wspace=.48,hspace=.48)
selected = []
for col,(species,pair,score) in enumerate([('Bacillus',('Factor2','Factor6'),13.700353),
                                         ('Rhodanobacter',('Factor2','Factor1'),12.988972)]):
    key=species.lower()
    r2=pd.read_csv(OUT/f'{key}_model_variance.csv')
    z=pd.read_csv(OUT/f'{key}_model_scores.csv')
    z['condition']=z.condition.replace({'Kana':'K','Kanamycin':'K','Phage':'P'})
    assert set(z.condition).issubset(set(old.COND_ORDER)), set(z.condition)
    w=pd.read_csv(OUT/f'{key}_model_weights.csv')
    ax=axes[0,col]
    vals=r2[old.OMICS].values
    im=ax.imshow(vals,aspect='auto',cmap='YlGnBu',vmin=0,vmax=50)
    sig=old.read_sig_factors(old.species_cfg[species]['anova'])
    ax.set_yticks(range(len(r2)),[f+('*' if f in sig else '') for f in r2.factor])
    ax.set_xticks(range(3),old.OMICS); ax.xaxis.tick_top()
    for i in range(10):
        for j in range(3):
            ax.text(j,i,f'{vals[i,j]:.1f}',ha='center',va='center',fontsize=14,
                    color='white' if vals[i,j]>26 else '#202020')
    ax.set_title(species+'\nVariance explained by factor and omics layer',loc='left',pad=25,fontsize=19)
    ax.text(0,-.1,'* Condition association: FDR < 0.05',transform=ax.transAxes,fontsize=14)
    ax=axes[1,col]
    for cond in old.COND_ORDER:
        sub=z[z.condition==cond]
        if sub.empty: continue
        xy=sub[list(pair)].values
        color=old.COND_COLORS[cond]
        ax.scatter(xy[:,0],xy[:,1],s=65,c=color,edgecolor='white',zorder=3)
        if len(xy)==3: ax.add_patch(Polygon(old.make_replicate_triangle(xy),color=color,alpha=.18))
        ax.annotate({'K':'Kan','P':'Ph'}.get(cond,cond),xy.mean(axis=0),xytext=(7,7),
                    textcoords='offset points',fontsize=14,fontweight='bold')
    ax.axhline(0,ls='--',lw=.7,c='#bac4d0'); ax.axvline(0,ls='--',lw=.7,c='#bac4d0')
    ax.margins(.18); ax.set_xlabel(pair[0]); ax.set_ylabel(pair[1])
    ax.set_title(f'{species}: {pair[0]} versus {pair[1]}\nPair-separation score = {score:.2f}',
                 loc='left',fontsize=19,pad=18)
    times='8 h; late = 24 h' if col==0 else '24 h; late = 48 h'
    ax.text(0,-.2,'Untreated sampling times: mid = '+times,transform=ax.transAxes,fontsize=14)
    z[['sample','condition',*pair]].to_csv(OUT/f'{key}_plotted_scores.csv',index=False)
    ax=axes[2,col]
    chunks=[]
    for view in old.OMICS:
        for factor in pair:
            sub=old.select_features(w,factor,view).copy()
            sub['display']=sub.feature.map(lambda f:old.format_feature(f,view,species))
            sub=sub[sub.display.notna()]
            sub['display']=sub['display'].map(lambda s:textwrap.fill(s.replace('\n',''),width=36))
            chunks.append(sub)
    chosen=pd.concat(chunks); chosen['species']=species; selected.append(chosen)
    rows=list(dict.fromkeys(zip(chosen.view,chosen.feature,chosen.display)))
    maxval=chosen.abs_weight.max()
    for i,(view,feature,label) in enumerate(rows):
        y=len(rows)-1-i
        ax.axhspan(y-.5,y+.5,color={'Transcriptome':'#f4f8ff','Proteome':'#fffaf4','Metabolome':'#f5fff1'}[view],zorder=0)
        for c,(factor,v) in enumerate([(f,v) for f in pair for v in old.OMICS]):
            sub=chosen[(chosen.factor==factor)&(chosen.view==v)&(chosen.feature==feature)]
            if view==v and len(sub):
                row=sub.iloc[0]
                ax.scatter(c,y,s=40+230*row.abs_weight/maxval,c=old.SIGN_COLORS[row['sign']],edgecolor='#555',lw=.4)
    ax.set_yticks(range(len(rows)),[x[2] for x in rows][::-1],fontsize=13)
    ax.set_xticks(range(6),[f'{f}\n{v}' for f in pair for v in old.OMICS],fontsize=13)
    ax.set_xlim(-.5,5.5); ax.set_ylim(-.5,len(rows)-.5)
    ax.axvline(2.5,c='#aeb8c4',lw=1)
    ax.set_title(species+' factor loadings\nBubble area scales with |loading|; red positive, blue negative',loc='left',fontsize=18,pad=20)
    ax.tick_params(length=0)
    for spine in ax.spines.values(): spine.set_visible(False)
for letter,ax in zip('ABCDEF',axes.flat):
    ax.text(-.17,1.12,letter,transform=ax.transAxes,fontsize=26,fontweight='bold')
cax=fig.add_axes([.965,.754,.012,.14]); fig.colorbar(im,cax=cax,label='Variance explained (%)')
pd.concat(selected).to_csv(OUT/'plotted_loadings.csv',index=False)
fig.savefig(OUT/'Figure4_verified.pdf',bbox_inches='tight')
fig.savefig(OUT/'Figure4_verified.svg',bbox_inches='tight')
fig.savefig(OUT/'Figure4_verified.png',dpi=400,bbox_inches='tight')
fig.savefig(OUT/'Figure4_preview.png',dpi=70,bbox_inches='tight')
