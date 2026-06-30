#!/usr/bin/env python
"""
CR vs non-CR TCRdist clustering within mutation-associated TCR sets,
CD8 and CD4 analysed SEPARATELY, two statistics:
  - mean pairwise CDR3 TCRdist
  - mean nearest-neighbour CDR3 TCRdist
For each candidate-neoantigen mutation with >=5 cancer-reactive (CR)
associated TCRs in a chain (CR defined per lineage from CD8_CAT_final.tsv /
CD4_CAT_final.tsv), split associated TCRs into CR vs non-CR and compare each
group's within-set distance to B=1000 length- & frequency-matched controls.
Alpha and beta separately.
"""
import os, numpy as np, pandas as pd
from tcrdist.repertoire import TCRrep
from scipy import stats

OUT="/Users/wsun/research/GitHub/neo-TCR/somatic_mutation/code_TCGA/step8b_NetMHCpan/step8h_tcrdist"
CATd="/Users/wsun/research/CAT/data"
B=1000; MIN_SET=5; MIN_CAND=15; SEED=2026
rng=np.random.default_rng(SEED)
VJ={"alpha":("v_a_gene","j_a_gene","cdr3_a_aa","TRAV1-1*01","TRAJ1*01"),
    "beta": ("v_b_gene","j_b_gene","cdr3_b_aa","TRBV2*01","TRBJ1-1*01")}

def cdr3_distmat(seqs, chain):
    vg,jg,cg,vd,jd=VJ[chain]
    df=pd.DataFrame({cg:seqs,vg:vd,jg:jd,"count":1})
    tr=TCRrep(cell_df=df,organism="human",chains=[chain],compute_distances=True,deduplicate=False)
    return getattr(tr,f"pw_cdr3_{chain[0]}_aa").astype(float)

def set_stats(D, idx):
    sub=D[np.ix_(idx,idx)]; n=len(idx)
    mean_pw=sub[np.triu_indices(n,1)].mean()
    mean_nn=(sub+np.eye(n)*1e9).min(axis=1).mean()
    return mean_pw, mean_nn

def freq_bins(counts):
    edges=np.unique(np.quantile(counts,np.linspace(0,1,11)))
    return np.digitize(counts,edges[1:-1],right=False)

def matched_null(D, focal_idx, exclude_idx, ln, fb, logc, label):
    obs_pw, obs_nn = set_stats(D, focal_idx)
    excl=set(exclude_idx); cand=[]
    for i in focal_idx:
        L,bb=ln[i],fb[i]; sel=None
        for wlen in (0,1,2):
            for wbin in (0,1,2,3):
                m=(np.abs(ln-L)<=wlen)&(np.abs(fb-bb)<=wbin)
                c=[j for j in np.where(m)[0] if j not in excl]
                if len(c)>=MIN_CAND: sel=c; break
            if sel is not None: break
        if sel is None:
            c=[j for j in np.where(np.abs(ln-L)<=2)[0] if j not in excl]
            sel=c if c else [j for j in range(D.shape[0]) if j not in excl]
        cand.append(np.array(sel))
    npw=np.empty(B); nnn=np.empty(B); dlog=[]; lmatch=[]
    for b in range(B):
        chosen,cset=[],set()
        for k in rng.permutation(len(focal_idx)):
            avail=cand[k][~np.isin(cand[k],list(cset))] if cset else cand[k]
            if len(avail)==0: avail=cand[k]
            j=int(rng.choice(avail)); chosen.append(j); cset.add(j)
            if b==0: dlog.append(abs(logc[j]-logc[focal_idx[k]])); lmatch.append(int(ln[j]==ln[focal_idx[k]]))
        npw[b],nnn[b]=set_stats(D,chosen)
    return dict(group=label,n=len(focal_idx),
        obs_mean=obs_pw,null_mean=npw.mean(),z_mean=(obs_pw-npw.mean())/npw.std(),
        p_mean=(1+np.sum(npw<=obs_pw))/(B+1),
        obs_nn=obs_nn,null_nn=nnn.mean(),z_nn=(obs_nn-nnn.mean())/nnn.std(),
        p_nn=(1+np.sum(nnn<=obs_nn))/(B+1),
        len_match_rate=np.mean(lmatch),mean_dlogcount=np.mean(dlog))

# ---- load ----
assoc=pd.read_csv(f"{OUT}/step8k_assoc_all.tsv",sep="\t")
cand=set(pd.read_csv(f"{OUT}/../step8g_neoantigen_list.tsv",sep="\t").key)
assoc=assoc[assoc.key.isin(cand)]
cd8=pd.read_csv(f"{CATd}/CD8_CAT_final.tsv",sep="\t"); cd4=pd.read_csv(f"{CATd}/CD4_CAT_final.tsv",sep="\t")
CR={("CD8","alpha"):set(cd8.TRA_cdr3.dropna()),("CD8","beta"):set(cd8.TRB_cdr3.dropna()),
    ("CD4","alpha"):set(cd4.TRA_cdr3.dropna()),("CD4","beta"):set(cd4.TRB_cdr3.dropna())}
pools={"alpha":pd.read_csv(f"{OUT}/control_pool_alpha.tsv",sep="\t").drop_duplicates("cdr3").reset_index(drop=True),
       "beta": pd.read_csv(f"{OUT}/control_pool_beta.tsv",sep="\t").drop_duplicates("cdr3").reset_index(drop=True)}

# precompute distance matrices per chain (shared across lineages)
chaindat={}
for chain in ("alpha","beta"):
    pool=pools[chain]
    assert not (set(assoc[assoc.chain==chain].cdr3)-set(pool.cdr3))
    seqs=pool.cdr3.to_numpy(); cnt=pool.carrier_count.to_numpy(); ln=pool.len.to_numpy()
    print(f"[{chain}] CDR3 TCRdist matrix {len(seqs)}x{len(seqs)} ...",flush=True)
    chaindat[chain]=dict(seqs=seqs,ln=ln,fb=freq_bins(cnt),logc=np.log(cnt),
                         idx={s:i for i,s in enumerate(seqs)},D=cdr3_distmat(seqs,chain))

rows=[]
for lineage in ("CD8","CD4"):
    for chain in ("alpha","beta"):
        cd=chaindat[chain]; crset=CR[(lineage,chain)]
        a=assoc[assoc.chain==chain]
        for key,g in a.groupby("key"):
            seqs_mut=g.cdr3.unique()
            cr=[cd["idx"][s] for s in seqs_mut if s in crset]
            nc=[cd["idx"][s] for s in seqs_mut if s not in crset]
            if len(cr)<MIN_SET: continue
            excl=cr+nc
            r=matched_null(cd["D"],cr,excl,cd["ln"],cd["fb"],cd["logc"],"CR"); r.update(lineage=lineage,key=key,chain=chain); rows.append(r)
            line=f"[{lineage} {chain:5s}] {key:22s} CR n={len(cr):2d} z_mean={r['z_mean']:+5.2f} z_nn={r['z_nn']:+5.2f}"
            if len(nc)>=MIN_SET:
                r2=matched_null(cd["D"],nc,excl,cd["ln"],cd["fb"],cd["logc"],"nonCR"); r2.update(lineage=lineage,key=key,chain=chain); rows.append(r2)
                line+=f" | nonCR n={len(nc):2d} z_mean={r2['z_mean']:+5.2f} z_nn={r2['z_nn']:+5.2f}"
            print(line,flush=True)

res=pd.DataFrame(rows)[["lineage","key","chain","group","n","obs_mean","null_mean","z_mean","p_mean",
    "obs_nn","null_nn","z_nn","p_nn","len_match_rate","mean_dlogcount"]]
res.to_csv(f"{OUT}/step8l_tcrdist_CR_results.tsv",sep="\t",index=False)

print("\n================ SUMMARY (z<0 = closer than matched control) ================")
for lineage in ("CD8","CD4"):
    for chain in ("alpha","beta"):
        d0=res[(res.lineage==lineage)&(res.chain==chain)]
        if d0.empty: print(f"[{lineage} {chain}] no mutations with >=5 CR"); continue
        for stat,zc in (("mean","z_mean"),("nn","z_nn")):
            for grp in ("CR","nonCR"):
                d=d0[d0.group==grp]; m=len(d)
                if m==0: continue
                closer=int((d[zc]<0).sum()); sp=stats.binomtest(closer,m,0.5,alternative="greater").pvalue
                print(f"[{lineage} {chain:5s} {stat:4s} {grp:5s}] {m:2d} muts | closer {closer}/{m} (sign p={sp:.2g}) | mean z={d[zc].mean():+.3f}")
            w=d0.pivot_table(index="key",columns="group",values=zc).dropna()
            if len(w)>=5 and "CR" in w and "nonCR" in w:
                ps=stats.wilcoxon(w["CR"],w["nonCR"]).pvalue
                print(f"   -> paired CR vs nonCR ({stat}, n={len(w)}): CR={w['CR'].mean():+.3f} nonCR={w['nonCR'].mean():+.3f} Wilcoxon p={ps:.3g}")
print("\nwrote step8l_tcrdist_CR_results.tsv")
