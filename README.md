# neo-TCR: neoantigen-specific TCRs

Identification and characterization of neoantigen-specific TCRs and the corresponding neoantigens

---

## Repository Layout

- **classify_TCR_by_NN**  
  Classify cancer-reactive TCRs using scRNA-/scTCR-seq data and neural networks trained on cell-level labels.

- **data_McPAS**  
  TCR annotation data from [McPAS-TCR](http://friedmanlab.weizmann.ac.il/McPAS-TCR/).

- **data_TCR**  
  Curated datasets from earlier studies reporting cancer-reactive TCRs (typically <100 per study), with R scripts for processing.

- **data_VDJdb**  
  TCR annotation data from [VDJdb](https://vdjdb.cdr3.net/).

- **somatic_mutation**  
  Analyses of associations between somatic mutations and TCR repertoires; predictors of mutation status built from TCR profiles.

---

## Filtering Cancer-Reactive T Cells

- [zz1_prepare_data.Rmd](https://github.com/Sun-lab/CAT/blob/c17bcc4b17bfc1c25196b22260a69c76df286149/zz1_prepare_data.Rmd) and [zz2_construct_training_data.Rmd](https://github.com/Sun-lab/CAT/blob/main/zz2_construct_training_data.Rmd): select cancer-reactive T cells based on signature scores, cluster, clone size, and antigen annotation.  
- `classify_cancer_reactive_T/step3_update_{CD4,CD8}_anno`: add neural network prediction results.  
- `somatic_mutation/code_TCGA/step6_cancer_reactive_T_signature.Rmd`: export signatures for mutation–TCR analysis.

---
