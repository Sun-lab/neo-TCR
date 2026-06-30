#!/usr/bin/env python
"""
TCRdist clustering test for candidate-neoantigen TCRs.

For each candidate neoantigen mutation (HLA-I, step8g), test whether its
associated TCRs (alpha / beta analysed separately, CDR3-only TCRdist) are
closer to each other than expected. The null is built from control TCR sets
drawn from the tested-TCR universe, matched to the associated set on both
population frequency (carrier count) and CDR3 length. Distances are the
CDR3 component of TCRdist (tcrdist3 default params: gap=4, ntrim=3, ctrim=2).
"""
import os, numpy as np, pandas as pd
from tcrdist.repertoire import TCRrep

OUT = "/Users/wsun/research/GitHub/neo-TCR/somatic_mutation/code_TCGA/step8b_NetMHCpan/step8h_tcrdist"
B        = 1000      # control draws per mutation
MIN_SET  = 5         # minimum associated TCRs per chain
MIN_CAND = 15        # minimum matched candidates per associated TCR (else widen)
SEED     = 2026
rng = np.random.default_rng(SEED)

VJ = {"alpha": ("v_a_gene","j_a_gene","cdr3_a_aa","TRAV1-1*01","TRAJ1*01"),
      "beta":  ("v_b_gene","j_b_gene","cdr3_b_aa","TRBV2*01",  "TRBJ1-1*01")}

def cdr3_distmat(seqs, chain):
    """CDR3-only TCRdist matrix among `seqs` (dummy constant V/J)."""
    vg, jg, cg, vd, jd = VJ[chain]
    df = pd.DataFrame({cg: seqs, vg: vd, jg: jd, "count": 1})
    tr = TCRrep(cell_df=df, organism="human", chains=[chain],
                compute_distances=True, deduplicate=False)
    return getattr(tr, f"pw_cdr3_{chain[0]}_aa").astype(float)

def set_stats(D, idx):
    """mean pairwise dist and mean nearest-neighbour dist over D[idx,idx]."""
    sub = D[np.ix_(idx, idx)]
    n = len(idx); iu = np.triu_indices(n, 1)
    mean_pw = sub[iu].mean()
    off = sub + np.eye(n)*1e9
    mean_nn = off.min(axis=1).mean()
    return mean_pw, mean_nn

def freq_bins(counts):
    """decile bins of carrier_count (population frequency numerator)."""
    edges = np.unique(np.quantile(counts, np.linspace(0,1,11)))
    return np.digitize(counts, edges[1:-1], right=False)

def run_chain(chain, assoc, pool):
    cg = VJ[chain][2]
    pool = pool.reset_index(drop=True)
    seqs = pool["cdr3"].to_numpy()
    cnt  = pool["carrier_count"].to_numpy()
    ln   = pool["len"].to_numpy()
    fb   = freq_bins(cnt)
    idx_of = {s:i for i,s in enumerate(seqs)}
    logc = np.log(cnt)

    print(f"[{chain}] computing {len(seqs)}x{len(seqs)} CDR3 TCRdist matrix ...", flush=True)
    D = cdr3_distmat(seqs, chain)

    rows = []
    muts = assoc[assoc.chain==chain].groupby("key_pro")["cdr3"].apply(lambda s: list(dict.fromkeys(s)))
    for kp, aset in muts.items():
        a_idx = [idx_of[s] for s in aset if s in idx_of]
        if len(a_idx) < MIN_SET:
            continue
        a_set = set(a_idx)
        obs_pw, obs_nn = set_stats(D, a_idx)

        # per-associated-TCR matched candidate pools (exclude this mut's own set)
        cand = []
        for i in a_idx:
            L, bb = ln[i], fb[i]
            sel = None
            for wlen in (0,1,2):
                for wbin in (0,1,2,3):
                    m = (np.abs(ln-L)<=wlen) & (np.abs(fb-bb)<=wbin)
                    c = [j for j in np.where(m)[0] if j not in a_set]
                    if len(c) >= MIN_CAND:
                        sel = c; break
                if sel is not None: break
            if sel is None:  # last resort: same length any freq, else any
                c = [j for j in np.where(np.abs(ln-L)<=2)[0] if j not in a_set]
                sel = c if c else [j for j in range(len(seqs)) if j not in a_set]
            cand.append(np.array(sel))

        # B matched control draws (distinct sequences within a draw)
        null_pw = np.empty(B); null_nn = np.empty(B)
        dlog = []; lmatch = []
        for b in range(B):
            order = rng.permutation(len(a_idx))
            chosen, cset = [], set()
            for k in order:
                avail = cand[k][~np.isin(cand[k], list(cset))] if len(cset) else cand[k]
                if len(avail)==0:
                    avail = cand[k]
                j = int(rng.choice(avail))
                chosen.append(j); cset.add(j)
                if b==0:
                    dlog.append(abs(logc[j]-logc[a_idx[k]])); lmatch.append(int(ln[j]==ln[a_idx[k]]))
            null_pw[b], null_nn[b] = set_stats(D, chosen)

        p_pw = (1+np.sum(null_pw <= obs_pw))/(B+1)
        p_nn = (1+np.sum(null_nn <= obs_nn))/(B+1)
        z_pw = (obs_pw-null_pw.mean())/null_pw.std()
        rows.append(dict(key_pro=kp, chain=chain, n_tcr=len(a_idx),
            obs_mean=obs_pw, null_mean=null_pw.mean(), null_sd=null_pw.std(),
            null_lo=np.quantile(null_pw,.025), null_hi=np.quantile(null_pw,.975),
            z=z_pw, p_emp=p_pw, obs_nn=obs_nn, null_nn_mean=null_nn.mean(), p_nn=p_nn,
            len_match_rate=np.mean(lmatch), mean_dlogcount=np.mean(dlog)))
        print(f"  {kp:20s} n={len(a_idx):3d} obs={obs_pw:6.2f} null={null_pw.mean():6.2f} "
              f"z={z_pw:+5.2f} p={p_pw:.3f}", flush=True)
    return pd.DataFrame(rows)

# ---- load ----
assoc = pd.read_csv(os.path.join(OUT,"assoc_tcrs.tsv"), sep="\t")
pool_a = pd.read_csv(os.path.join(OUT,"control_pool_alpha.tsv"), sep="\t")
pool_b = pd.read_csv(os.path.join(OUT,"control_pool_beta.tsv"),  sep="\t")

res = pd.concat([run_chain("alpha", assoc, pool_a),
                 run_chain("beta",  assoc, pool_b)], ignore_index=True)
res.to_csv(os.path.join(OUT,"step8i_tcrdist_results.tsv"), sep="\t", index=False)

# ---- aggregate ----
from scipy import stats
print("\n================ SUMMARY ================")
for ch in ("alpha","beta"):
    d = res[res.chain==ch]
    m = len(d); nsig = int((d.p_emp<0.05).sum()); ncloser = int((d.z<0).sum())
    sign_p = stats.binomtest(ncloser, m, 0.5, alternative="greater").pvalue
    wilc = stats.wilcoxon(d.obs_mean - d.null_mean, alternative="less").pvalue if m>5 else np.nan
    stouffer = d.z.sum()/np.sqrt(m)
    print(f"[{ch}] {m} mutations | closer-than-null (z<0): {ncloser}/{m} (sign p={sign_p:.3g}) | "
          f"emp p<0.05: {nsig} | mean z={d.z.mean():+.3f} Stouffer Z={stouffer:+.2f} | "
          f"Wilcoxon(obs<null) p={wilc:.3g}")
    print(f"      matching: len-match={d.len_match_rate.mean()*100:.0f}%  "
          f"mean|Dlog(count)|={d.mean_dlogcount.mean():.3f}")
print("\nwrote step8i_tcrdist_results.tsv")
