# TRACERx Somatic Mutation Prediction Model Testing

This directory contains code and workflows for **testing somatic mutation prediction models trained on TCGA data** using **TRACERx** datasets.

## Directory Overview

- **`somatic_mutation/`**  
  Build a **Sample × Somatic Mutation** matrix using mutation data from the TRACERx paper.

- **`rna-seq/`**  
  - Reconstruct TCRs from TRACERx RNA-seq data using **TRUST4**.  
  - Apply TCGA-trained models for evaluation.

- **`tcr-seq/`**  
  - Call TCR sequences from TRACERx TCR-seq data.  
  - Apply TCGA-trained models for evaluation.

- **`Compare_TCR_datasets/`**  
  Compare TCR repertoires obtained from **RNA-seq** and **TCR-seq** datasets.
