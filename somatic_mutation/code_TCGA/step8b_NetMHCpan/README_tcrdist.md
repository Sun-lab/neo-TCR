# TCRdist analysis of mutation-associated TCRs (steps 8g–8m)

Pipeline that (1) builds a candidate-neoantigen list, then asks whether the TCRs
associated with a somatic mutation are **more similar to each other** (smaller CDR3
TCRdist) than length- & frequency-matched control TCRs — overall, and split by
cancer-reactive (CR) vs non-CR status (CD8 / CD4 separately).

- **Distance metric:** CDR3-only TCRdist (tcrdist3 default params: gap=4, ntrim=3,
  ctrim=2). V-genes are unavailable in our data, so only the CDR3 component is used;
  alpha and beta are analysed separately (unpaired bulk TCRs).
- **Control / null:** for each TCR set, B=1000 control sets drawn from the eligible
  tested-TCR universe, matched per-TCR on CDR3 length and carrier frequency.
- **Associated TCRs:** full-TCGA predictors at pA=0.1 / pB=0.2
  (`TCGA_res_TCR_predictors_pval_pA_1.0e-01_pB_2.0e-01.rds`).

All code lives in this folder (`step8b_NetMHCpan/`); all data/figure outputs go to
`step8b_NetMHCpan/step8h_tcrdist/` (except the step8g lists, which sit one level up
here). Each `.R` writes a `.Rout` log beside it. **Python scripts need the tcrdist3
venv** at `…/scratchpad/tcrdist_env` (tcrdist3 + parasail; parasail needs autotools).

---

## Run order

| # | script | lang | reads | writes |
|---|--------|------|-------|--------|
| 8g | `step8g_neoantigen_list.R` | R | step8d carrier-freq table, TCR predictor RDS, mut-asso TSVs | `step8g_neoantigen_list.tsv`, `step8g_neoantigen_list_TCR.tsv` (in this folder) |
| 8h | `step8h_export_for_tcrdist.R` | R | step8g list, TCR predictor RDS, `TCR_{alpha,beta}_mut_asso_rd_x_with_ct.tsv` | `step8h_tcrdist/assoc_tcrs.tsv`, `control_pool_alpha.tsv`, `control_pool_beta.tsv` |
| 8i | `step8i_tcrdist_analysis.py` | py | `assoc_tcrs.tsv`, control pools | `step8i_tcrdist_results.tsv` |
| 8j | `step8j_tcrdist_plot.py` | py | `step8i_tcrdist_results.tsv` | `step8j_tcrdist_clustering.pdf` |
| 8k | `step8k_export_CR.R` | R | TCR predictor RDS, `CD8/CD4_CAT_final.tsv`, step8g list | `step8h_tcrdist/step8k_assoc_all.tsv`, `step8k_overlap_counts.tsv` |
| 8l | `step8l_tcrdist_CR.py` | py | `step8k_assoc_all.tsv`, `CD8/CD4_CAT_final.tsv`, control pools | `step8l_tcrdist_CR_results.tsv` |
| 8m | `step8m_plot_CR.py` | py | `step8l_tcrdist_CR_results.tsv` | `step8m_tcrdist_CR_vs_nonCR.pdf` |

---

## What each script does

- **8g — candidate neoantigen list.** From the step8d HLA-presentation table, keep
  (mutation, HLA) pairs where the allele presents ≥1 neo-peptide AND its normalized
  frequency in carriers (E) > 0.05. Adds associated TCRs + p-values in the `_TCR` file.
  → 58 pairs, **25 candidate mutations**.

- **8h — export for the basic clustering test.** For the 25 candidate mutations, the
  associated TCRs per chain (`assoc_tcrs.tsv`), plus the eligible control universe with
  each TCR's carrier count and CDR3 length (`control_pool_{alpha,beta}.tsv`).

- **8i — basic clustering test.** Per mutation/chain, mean pairwise CDR3 TCRdist of the
  associated TCRs vs B=1000 matched-control sets (z, empirical p). Result: associated
  TCRs are **not** closer than matched controls.

- **8j — plot of 8i.** Per-mutation observed mean distance vs the matched-null 95% band,
  α/β panels. (No CD8/CD4 dimension — TCGA-associated TCRs are not lineage-typed.)

- **8k — CR intersection + export.** CR TCR sets = unique TRA/TRB CDR3 from
  `CD8_CAT_final.tsv` / `CD4_CAT_final.tsv`. Exports associated TCRs for ALL mutations
  with an `is_CR` flag, plus per-mutation/chain CR-overlap counts. (Diagnostics: among
  the 25 candidates, ≥5 CR-associated → α 23 / β 8 with the CD8∪CD4 union.)

- **8l — CR vs non-CR clustering, CD8/CD4 separately.** For each candidate mutation with
  ≥5 CR-associated TCRs in a chain (CR defined per lineage), split associated TCRs into
  CR vs non-CR and compare each to matched controls, with **two statistics: mean pairwise
  and nearest-neighbour** TCRdist. Qualifying: CD8 α=22 β=7, CD4 α=16, CD4 β=0.
  Result: CR groups = matched null; non-CR α groups are more dispersed than matched
  controls (CD8 α mean p≈0.01).

- **8m — plot of 8l.** Paired CR vs non-CR z per mutation, 2 statistics × 3 analyzable
  lineage/chain panels.

---

## Output files (`step8h_tcrdist/`)

| file | description |
|------|-------------|
| `assoc_tcrs.tsv` | associated TCRs for the 25 candidate mutations (`key, key_pro, chain, cdr3`) |
| `control_pool_alpha.tsv`, `control_pool_beta.tsv` | eligible control TCRs (`cdr3, carrier_count, len`) |
| `step8i_tcrdist_results.tsv` | per-mutation/chain obs/null/z/p (basic test) |
| `step8j_tcrdist_clustering.pdf` | basic-test figure |
| `step8k_assoc_all.tsv` | associated TCRs, all 4,610 mutations (`key, chain, cdr3, is_CR`) |
| `step8k_overlap_counts.tsv` | per-mutation/chain `n_assoc, n_CR, n_nonCR` |
| `step8l_tcrdist_CR_results.tsv` | per-mutation `lineage, chain, group(CR/nonCR), z_mean, p_mean, z_nn, p_nn`, … |
| `step8m_tcrdist_CR_vs_nonCR.pdf` | CR-vs-non-CR figure |

---

## Known overlaps (intentional; not refactored)

- `cdr3_distmat()`, `set_stats()`, `freq_bins()` and the matched-null sampling logic are
  duplicated between `step8i` and `step8l` (`step8l` is the generalized version: adds the
  nearest-neighbour statistic, lineage loop, and CR/non-CR groups).
- `assoc_tcrs.tsv` (8h, 25 mutations, + `key_pro` label) is a subset of
  `step8k_assoc_all.tsv` (8k, all mutations, + `is_CR`). Both are kept because `step8i`
  reads the former and `step8h` also produces the control pools used downstream.

## Caveats

CDR3-only TCRdist (no V-gene); α/β separate (unpaired bulk); CR set from CD8/CD4
`_CAT_final.tsv` (= `label==1 & Prediction>0.5`); p-values are nominal (no multiple-test
correction across the per-group panels).
