# neo-TCR

**Characterization of Cancer-Reactive T Cells and Neoantigen-specific T Cell Receptors**

Wei Sun, Si Liu, Bo Yu, Zhehao Zhang, Philip Bradley, Marie Bleakley

This repository contains the code and curated data for identifying and characterizing
**cancer-reactive (CR) T cells** and the **neoantigen-specific T cell receptors (TCRs)** and
**neoantigens** they recognize.

> Companion repository: the data-preparation steps for the cancer-reactive T cell atlas live in
> [**Sun-lab/CAT**](https://github.com/Sun-lab/CAT). Together, `CAT` and `neo-TCR` implement the two
> computational platforms described in the paper.

---

## Overview

Not all tumor-infiltrating T cells are cancer-reactive; many are bystanders. This project develops
two complementary, sequencing-based computational platforms:

### 1. CAT — Cancer-Associated T cells
Harmonizes cancer-reactive gene expression signatures from prior studies and applies them to an
**atlas of ~1 million T cells** (≈500k CD8⁺ and ≈500k CD4⁺) from **383 patients across 11 cancer
types**, pooled from five single-cell studies. A four-step workflow identifies high-confidence CR
T cells:

1. **Initial CR set** — cells with high positive and low negative signature scores.
2. **Intersect** with study-specific CR clusters (clusters enriched for the above).
3. **Refine** — keep cells from clonally expanded TCRs and drop TCRs annotated to non-human
   antigens (via VDJdb).
4. **Neural-network filter** — an autoencoder-based classifier scores genome-wide expression;
   cells with score > 0.5 are retained.

This yields **62,188 CD8⁺ and 7,322 CD4⁺ CR T cells** (~70k total) with per-cell signature scores
and TCR sequences. The baseline abundance of CD8⁺ CR T cells predicts response to immune
checkpoint inhibitor therapy.

### 2. Neo-TCR — neoantigen-specific TCRs and neoantigens
A high-throughput framework that screens for **public neoantigens** (shared across patients) and
their cognate TCRs by association analysis. Using TCR and somatic mutation data for **8,463 TCGA
patients**, each TCR–mutation pair is tested with logistic regression. Somatic mutations whose
status can be predicted from their associated TCRs are more likely to encode public neoantigens;
candidate neoantigen peptides are then prioritized by predicted HLA presentation. Predictions are
validated by cross-validation and replicated in two independent datasets (Zheng et al. 2021 and
TRACERx). This prioritizes **47 candidate neoantigen-encoding mutations**, 8 of which have prior
experimental support.

A key convergence result: T cells carrying mutation-associated (Neo-TCR) TCRs are enriched for CR
gene expression signatures — integrating the two platforms improves T cell–based biomarkers for
cancer detection and monitoring.

---

## Repository layout

| Path | Platform | Description |
| --- | --- | --- |
| [`classify_cancer_reactive_T/`](classify_cancer_reactive_T) | CAT | End-to-end workflow that labels CR vs. non-CR T cells from signatures and a neural-network classifier. |
| [`classify_TCR_by_NN/`](classify_TCR_by_NN) | CAT ↔ Neo-TCR | Tests whether CR status can be predicted from **TCR sequence alone** (CNN ensemble); motivates the association-based Neo-TCR approach. |
| [`somatic_mutation/`](somatic_mutation) | Neo-TCR | TCR–somatic-mutation association analysis, mutation-prediction models, validation, and neoantigen prioritization. |
| [`data_TCR/`](data_TCR) | reference | Experimentally identified cancer-reactive TCRs curated from prior studies. |
| [`data_McPAS/`](data_McPAS) | reference | TCR–antigen annotations from [McPAS-TCR](http://friedmanlab.weizmann.ac.il/McPAS-TCR/). |
| [`data_VDJdb/`](data_VDJdb) | reference | TCR–antigen annotations from [VDJdb](https://vdjdb.cdr3.net/) (used to drop non-human-antigen TCRs). |
| [`conda_setup.sh`](conda_setup.sh) | — | Creates the `careTCR` conda environment used by the Python notebooks. |

---

## CAT: identifying cancer-reactive T cells

### `classify_cancer_reactive_T/`
The main CAT workflow (Python / Jupyter notebooks):

| Step | File(s) | Purpose |
| --- | --- | --- |
| 1 | `step1a_prepare_CD8.ipynb`, `step1b_prepare_CD4.ipynb` | Assemble and QC the CD8⁺ / CD4⁺ T cell atlas and score gene signatures. |
| 2 | `step2_classify_careT.ipynb`, `step2_run_classify_careT.py` | Initial signature-based CR labeling and study-specific CR clustering. |
| NN | `CD8_<study>_lr_0.001_bs_32.ipynb`, `CD4_<study>_lr_0.001_bs_32.ipynb` | Train the autoencoder-based CR classifier per study (leave-one-study-out), for studies `Zheng_2021`, `Liu_2022`, `Liu_2025`, `Chow_2023`, `Chen_2024`. |
| 3 | `step3_update_CD8_anno.ipynb`, `step3_update_CD4_anno.ipynb` | Merge neural-network prediction scores into the cell annotation. |
| 4 | `step4_classify_careT_CD8.ipynb`, `step4_classify_careT_CD4.ipynb` | Produce the final high-confidence CR T cell sets (NN score > 0.5). |

Outputs: trained models (`model/`), per-study predictions (`classification_result/`), differential
expression (`DE_results/`), and figures (`figures/`).

### `classify_TCR_by_NN/`
Asks the complementary question: *can a TCR's amino-acid sequence alone tell us whether it came
from a CR T cell?* TCRs are labeled by the CR status of their cell of origin, and a CNN ensemble is
trained across a grid of inputs (α / β / both chains × CDR3 ± V gene), encodings (one-hot, BLOSUM62,
Atchley, PCA), CNN depths, and dense-layer sizes. Test AUC is averaged over a 20-model ensemble.
The modest accuracy of these models — TCR sequence alone does not reliably distinguish CR TCRs —
motivates the association-based strategy of **Neo-TCR**. See [`classify_TCR_by_NN/README.md`](classify_TCR_by_NN/README.md)
for the full feature/architecture grid and performance plots.

---

## Neo-TCR: associating TCRs with somatic mutations and neoantigens

### `somatic_mutation/code_TCGA/`
The core association and prediction pipeline (R / R Markdown). See
[`somatic_mutation/code_TCGA/README.md`](somatic_mutation/code_TCGA/README.md) for details.

| Step | Purpose |
| --- | --- |
| `step1_prepare_data` | Prepare TCGA TCR and somatic-mutation tables; keep TCRs/mutations seen in ≥5 patients; add cancer type. |
| `step2b_TCR_mut_asso_rd_interaction` | Test each TCR–mutation association by logistic regression with a likelihood-ratio test (modeling TCR read depth, mutation status, and their interaction); `step2_simulation` calibrates the test. |
| `step3_check_results`, `step3b_cancer_type` | Inspect associations and confirm cancer type is not a confounder. |
| `step4_build_prediction_model` | Aggregate associated TCRs into per-mutation predictors of mutation status. |
| `cross_validation/` | Leave-one-cancer-type-out cross-validation; AUC/FDR estimation and p-value-threshold selection. |
| `step5_check_Zheng_et_al` | External validation on Zheng et al. 2021 scTCR data (e.g., BRAF V600E). |
| `step6_cancer_reactive_T_signature` | Show that cells with mutation-associated TCRs carry CR signatures (links Neo-TCR back to CAT). |
| `step7_extract_peptides` | Extract candidate neoantigen peptides around predictable mutations and characterize them. |
| `step8_check_HLA_overlap`, `step8b_NetMHCpan` | Prioritize peptides by predicted HLA presentation (NetMHCpan). |
| `step9_check_known_neoantigens` | Cross-reference candidates against experimentally validated neoantigens. |

### `somatic_mutation/code_TRACERx/`
Independent replication of the TCGA models on the lung **TRACERx** cohort:
- `somatic_mutation/` — build the sample × mutation matrix from TRACERx calls.
- `rna-seq/` — reconstruct TCRs from RNA-seq with **TRUST4**, then apply TCGA-trained models.
- `tcr-seq/` — call TCRs from TCR-seq with **Decombinator**, then apply TCGA-trained models.
- `Compare_TCR_datasets.Rmd` — compare RNA-seq vs. TCR-seq repertoires.

### `somatic_mutation/tcr_cluster_function/`
A reusable function, `get_predictor_cluster()`, that scores a new individual's CDR3 repertoire
against the mutation-associated TCR clusters (for `alpha`/`beta` chains at height cutoffs 12 or 20).
See [`somatic_mutation/tcr_cluster_function/README.md`](somatic_mutation/tcr_cluster_function/README.md)
for inputs, outputs, and a runnable example.

---

## Reference / annotation data

- **`data_TCR/`** — cancer-reactive TCRs isolated in experimental studies (Caushi 2021, Hanada 2022,
  Lowery 2022, Sudmeier 2022, Veatch 2022); `step1_prepare_TCR.Rmd` harmonizes them.
- **`data_McPAS/`** — McPAS-TCR annotations (downloaded 2023-09-25).
- **`data_VDJdb/`** — VDJdb annotations (releases 2023-06-01 and 2025-02-21).

---

## Data availability

Several inputs are access-controlled and are **not** redistributed here:

- **TCGA TCRs** — `TCGA_mitcr_cdr3_result_161008.tsv` from the
  [PanImmune](https://gdc.cancer.gov/about-data/publications/panimmune) study.
- **TCGA somatic mutations** — `mc3.v0.2.8.CONTROLLED.maf.gz` from the
  [PanCanAtlas](https://gdc.cancer.gov/about-data/publications/pancanatlas).
- **Single-cell RNA-seq / TCR-seq** — GEO accessions
  [GSE156728](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156728) (pan-cancer, Zheng 2021),
  [GSE179994](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179994) (lung, Liu 2022),
  [GSE212217](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212217) (endometrial, Chow 2023),
  [GSE236581](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236581) (colorectal, Chen 2024),
  [GSE243013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243013) (lung, Liu 2025).
- **TRACERx** — TCR-seq (SRA), RNA-seq ([EGAD00001009862](https://ega-archive.org/datasets/EGAD00001009862)),
  and published somatic mutations.

---

## Setup

The Python notebooks (CAT and `classify_TCR_by_NN`) use a conda environment:

```bash
bash conda_setup.sh        # creates and configures the `careTCR` env
conda activate careTCR
```

The environment includes `numpy`, `pandas`, `matplotlib`, `seaborn`, `scikit-learn`, `tensorflow`,
and `notebook` (Python 3.11). The `somatic_mutation/` and `data_*/` analyses are written in R / R
Markdown; external tools used in the TRACERx validation include **TRUST4**, **Decombinator**, and
**NetMHCpan**.

---

## Citation

If you use this code or data, please cite:

> Sun W, Liu S, Yu B, Zhang Z, Bradley P, Bleakley M.
> *Characterization of Cancer-Reactive T Cells and Neoantigen-specific T Cell Receptors.*

Related software: [Sun-lab/CAT](https://github.com/Sun-lab/CAT) ·
[Sun-lab/neo-TCR](https://github.com/Sun-lab/neo-TCR)

## License

See [`LICENSE`](LICENSE) (Academic and Research Use License).
