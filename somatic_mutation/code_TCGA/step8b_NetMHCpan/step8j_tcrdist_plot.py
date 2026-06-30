#!/usr/bin/env python
"""Plot observed within-set CDR3 TCRdist vs frequency/length-matched null band."""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT="/Users/wsun/research/GitHub/neo-TCR/somatic_mutation/code_TCGA/step8b_NetMHCpan/step8h_tcrdist"
r=pd.read_csv(os.path.join(OUT,"step8i_tcrdist_results.tsv"),sep="\t")

fig,axes=plt.subplots(1,2,figsize=(9.5,5.4))
for ax,ch in zip(axes,("alpha","beta")):
    d=r[r.chain==ch].sort_values("obs_mean").reset_index(drop=True)
    y=np.arange(len(d))
    # matched-null 95% band + mean
    for i,row in d.iterrows():
        ax.plot([row.null_lo,row.null_hi],[i,i],color="0.78",lw=4,solid_capstyle="round",zorder=1)
        ax.plot(row.null_mean,i,marker="|",color="0.45",ms=9,zorder=2)
    # observed: red=below band (closer), blue=above (farther), dark=inside
    col=np.where(d.obs_mean<d.null_lo,"#d6604d",np.where(d.obs_mean>d.null_hi,"#4393c3","#222222"))
    ax.scatter(d.obs_mean,y,c=col,s=34,zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(d.key_pro,fontsize=7)
    ax.set_xlabel("CDR3 TCRdist (mean pairwise distance)",fontsize=10)
    ncloser=int((d.z<0).sum()); m=len(d)
    ax.set_title(f"{ch}  ({ncloser}/{m} closer than null; mean z = {d.z.mean():+.2f})",fontsize=10)
    ax.margins(y=0.02)
    ax.grid(axis="x",color="0.92")
# shared legend below the panels (kept off the data area)
from matplotlib.lines import Line2D
leg=[Line2D([0],[0],marker="o",color="w",markerfacecolor="#222222",ms=7,label="observed (inside band)"),
     Line2D([0],[0],marker="o",color="w",markerfacecolor="#d6604d",ms=7,label="observed below band (closer)"),
     Line2D([0],[0],marker="o",color="w",markerfacecolor="#4393c3",ms=7,label="observed above band (farther)"),
     Line2D([0],[0],color="0.78",lw=4,label="matched-null 95% band"),
     Line2D([0],[0],marker="|",color="0.45",ms=9,lw=0,label="null mean")]
fig.legend(handles=leg,fontsize=8.5,loc="lower center",ncol=3,frameon=False,
           bbox_to_anchor=(0.5,0.045))
fig.suptitle("Associated TCRs are not closer to each other than frequency- & length-matched controls",
             fontsize=11,y=0.995)
fig.text(0.5,0.012,"lower distance = TCRs more similar to each other",ha="center",fontsize=8,color="0.4")
fig.tight_layout(rect=[0,0.13,1,0.975])
f=os.path.join(OUT,"step8j_tcrdist_clustering.pdf")
fig.savefig(f); print("wrote",f)
