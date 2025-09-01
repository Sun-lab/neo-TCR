# TRACERx RNA-seq TCR Reconstruction & TCGA Model Validation

This workflow reconstructs TCRs from TRACERx RNA-seq data using **TRUST4** and validates somatic mutation prediction models trained on TCGA data.

## Data Source
- **RNA-seq data:** [EGAD00001009862](https://ega-archive.org/datasets/EGAD00001009862)
- **TCR reconstruction tool:** TRUST4

## Workflow Steps

### 01_run_trust4
- **Input:** TRACERx RNA-seq FASTQ files  
- **Output:** TRUST4 report of reconstructed TCRs

### 02_get_tcr
- **Input:** TRUST4 report of reconstructed TCRs  
- **Output:** Tumor (Patient) × TCR matrix

### 03_validate_TCGA_results
- **Input:** TRACERx somatic mutations & TCR data, TCGA prediction models  
- **Output:** Predictive AUC of TCGA models on TRACERx data
