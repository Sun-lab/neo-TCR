# Assessing Associations Between Somatic Mutations and TCR Using TCGA Data

This repository contains **R code** for analyzing the associations between somatic mutations and T cell receptors (TCR) using TCGA data.

> **Note:**  
> The dataset was obtained from [A Pan-Cancer Immunogenomic Analysis from TCGA](https://gdc.cancer.gov/about-data/publications/panimmune).  
> Due to access restrictions, we cannot share the data in this public repository.

---

## Workflow Overview

### 1. `step1_prepare_data` — Data Preparation
- Prepare TCR and somatic mutation datasets.
- Observed: Total TCR count per sample is strongly associated with T cell proportions.
- Select:
  - **2,333 TCR α chains**
  - **8,750 TCR β chains**  
  that appear in at least 5 subjects (out of 8,400+ subjects).
- Filter somatic mutations.
- Add cancer type information.
- Remove subjects with more than one cancer type.

---

### 2. `step2_TCR_mut_association` — Association Analysis
- Assess associations between TCRs and somatic mutations using **logistic regression**.
- Perform **Likelihood Ratio Tests (LRT)** for statistical significance.

---

### 3a. `step3a_check_results` — Results Review
- Inspect association results.
- Evaluate additional metrics to quantify TCR–mutation relationships.

---

### 3b. `step3b_cancer_type` — Cancer Type Effect
- Check for potential confounding effects of cancer type.
- **Conclusion:** No strong confounder effect from cancer type.

---

### 4. `step4_build_prediction_model` — Predictive Modeling
- Aggregate TCRs to build models predicting the presence of each somatic mutation.

---

### 5. `step5_check_Zheng_et_al` — External Validation
- Validate predictions using **Zheng et al.** dataset.

---

### 6. `step6_cancer_reactive_T_signature` — Cancer-Reactive Signature Analysis
- Assess cancer-reactive TCR signatures for TCRs significantly associated with somatic mutations.

---

### 7. `step7_extract_peptides` 
- Extract peptides around mutations for neoantigen analysis
- Charaterize the somatic mutations that are associted with TCRs

