all_4619_tested_mutations.tsv
  One row per tested mutation, with genomic/protein annotation, total carriers
  and eligibility for the current fixed-setting validation.

all_4619_mutations_in_sample_AUC_by_setting.tsv
  13,857 rows: all 4,619 mutations at each of three original alpha/beta TCR
  selection p-value threshold pairs (0.01/0.05, 0.05/0.10, 0.10/0.20).
  AUC and TCR-count columns retain the original exported values. Mutations
  absent from an original AUC export are retained, with AUC fields NA and an
  explicit source_AUC_record_available flag. Missing AUCs are not imputed.
  Original chain-specific n_alpha/n_beta can also be NA; values are preserved.


