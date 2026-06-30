#!/usr/bin/env python
"""Paired CR vs non-CR within-set TCRdist, CD8/CD4 separately, two statistics."""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

OUT="/Users/wsun/research/GitHub/neo-TCR/somatic_mutation/code_TCGA/step8b_NetMHCpan/step8h_tcrdist"
r=pd.read_csv(f"{OUT}/step8l_tcrdist_CR_results.tsv",sep="\t")

cols=[("CD8","alpha"),("CD8","beta"),("CD4","alpha")]   # CD4-beta has 0 qualifying mutations
rows=[("mean pairwise","z_mean"),("nearest-neighbour","z_nn")]
fig,axes=plt.subplots(len(rows),len(cols),figsize=(10,6.4),sharey="row")
for ri,(sname,zc) in enumerate(rows):
    for ci,(lin,ch) in enumerate(cols):
        ax=axes[ri,ci]
        d=r[(r.lineage==lin)&(r.chain==ch)]
        w=d.pivot_table(index="key",columns="group",values=zc).dropna()
        for _,row in w.iterrows():
            closer=row["CR"]<row["nonCR"]
            ax.plot([0,1],[row["CR"],row["nonCR"]],"-",color=("#4575b4" if closer else "#d73027"),alpha=0.5,lw=1,zorder=1)
            ax.scatter([0,1],[row["CR"],row["nonCR"]],s=20,color=["#2c7fb8","#999999"],zorder=2)
        for x,grp in zip((0,1),("CR","nonCR")):
            if grp in w: ax.plot([x-0.18,x+0.18],[w[grp].mean()]*2,color="black",lw=2.5,zorder=3)
        ax.axhline(0,ls="--",color="0.6",lw=1)
        p=stats.wilcoxon(w["CR"],w["nonCR"]).pvalue if (len(w)>=5 and "CR" in w and "nonCR" in w) else float("nan")
        ax.set_xticks([0,1]); ax.set_xticklabels(["CR","non-CR"])
        ax.set_xlim(-0.5,1.5); ax.grid(axis="y",color="0.93")
        crm=w["CR"].mean() if "CR" in w else float("nan"); ncm=w["nonCR"].mean() if "nonCR" in w else float("nan")
        ax.set_title(f"{lin} {ch} — {sname}\nn={len(w)}, CR z̄={crm:+.2f} / nonCR z̄={ncm:+.2f}, p={p:.2g}",fontsize=8.5)
        if ci==0: ax.set_ylabel(f"z ({sname}\nvs matched null)",fontsize=8.5)
fig.suptitle("CR vs non-CR mutation-associated TCR clustering (CD8 / CD4 separately)",fontsize=11)
fig.text(0.5,0.01,"lower z = TCRs more similar to each other than length/freq-matched controls; "
         "blue line = CR below non-CR",ha="center",fontsize=8,color="0.4")
fig.tight_layout(rect=[0,0.025,1,0.96])
f=f"{OUT}/step8m_tcrdist_CR_vs_nonCR.pdf"; fig.savefig(f); print("wrote",f)
