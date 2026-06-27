# =====================================================================
# Reviewer 3.6 (focused): hot vs cold performance on the PRIORITIZED
# mutations (the cv_optimal_pcut combos that pass the FDR filter and are
# reported in Fig 2C). Per-combo strata are tiny, so we POOL patients
# across prioritized mutations (same idea as s3's get_auc aggregator):
# one hot AUC and one cold AUC per chain, plus a depth-matched hot AUC.
# =====================================================================
suppressMessages({ library(data.table); library(pROC) })
source("../_utilities.R")

data_dir = "~/research/data/TCR/TCGA/_processed_data"
res_dir  = "~/research/data/TCR/TCGA/_results/_cross_validation"
set.seed(2026)
N_RAREFY = 50

sample_info = fread(file.path(data_dir, "sample_info_intersection.tsv"))
sample_info = sample_info[!is.na(T_cell_prop), ]
sample_info[, med_tcp := median(T_cell_prop), by = cancer_type]
sample_info[, hot_cold := ifelse(T_cell_prop > med_tcp, "hot", "cold")]

safe_auc <- function(y, x){
  if(length(unique(y)) < 2) return(NA_real_)
  as.numeric(auc(response = y, predictor = x, quiet = TRUE, direction = "<"))
}

pool_chain <- function(chain){
  combos = fread(file.path(res_dir,
            sprintf("cv_optimal_pcut_%s_mut_cancer_combos.tsv", chain)))
  TCR_data = fread(file.path(data_dir, sprintf("TCR_%s_5.tsv", chain)))
  mut_data = fread(file.path(data_dir, "mutation_freq_geq_5.tsv"))

  acc = data.table(label = integer(0), pred = numeric(0),
                   stratum = character(0), depth = numeric(0),
                   patient = character(0), cell = character(0))

  for(k in 1:nrow(combos)){
    c_type = combos$cancer[k]; mj = combos$mut[k]; pc = combos$pcut[k]
    r1 = fread(file.path(res_dir, sprintf("train_%s_wo_%s.tsv", chain, c_type)))

    si_c = sample_info[cancer_type == c_type, ]
    d1   = get_dat_one_cancer(c_type, TCR_data, mut_data, si_c)
    mut1 = d1$mut_mat; tcr1 = d1$TCR_mat
    if(!(mj %in% rownames(mut1))) next

    tcr_sel = r1[mut == mj & pval_LRT < pc, TCR]
    idx     = which(rownames(tcr1) %in% tcr_sel)
    if(length(idx) == 0) next
    pred = colSums(tcr1[idx, , drop = FALSE] > 0)
    y    = mut1[mj, ]

    sp = si_c[match(colnames(mut1), participant), ]
    acc = rbind(acc, data.table(label = as.integer(y), pred = as.numeric(pred),
          stratum = sp$hot_cold, depth = colSums(tcr1),
          patient = colnames(mut1), cell = paste(c_type, mj)))
  }
  acc[, chain := chain]
  acc
}

# Combined alpha+beta predictor (the "both" prioritized set, Fig 2C/2B).
# Mirrors s4_cross_validation_res_both.Rmd: concatenate alpha-TCRs (prefixed
# "alpha_") selected at pcut_A with beta-TCRs (prefixed "beta_") selected at
# pcut_B, then count detected clones. Patients pooled across combos.
pool_both <- function(){
  combos   = fread(file.path(res_dir, "cv_optimal_pcut_both_mut_cancer_combos.tsv"))
  TCR_a    = fread(file.path(data_dir, "TCR_alpha_5.tsv"))
  TCR_b    = fread(file.path(data_dir, "TCR_beta_5.tsv"))
  mut_data = fread(file.path(data_dir, "mutation_freq_geq_5.tsv"))

  acc = data.table(label = integer(0), pred = numeric(0),
                   stratum = character(0), depth = numeric(0),
                   patient = character(0), cell = character(0))

  for(k in 1:nrow(combos)){
    c_type = combos$cancer[k]; mj = combos$mut[k]
    pcA = combos$pcut_A[k];    pcB = combos$pcut_B[k]
    r_a = fread(file.path(res_dir, sprintf("train_alpha_wo_%s.tsv", c_type)))
    r_b = fread(file.path(res_dir, sprintf("train_beta_wo_%s.tsv",  c_type)))
    r_a[, TCR := paste0("alpha_", TCR)]
    r_b[, TCR := paste0("beta_",  TCR)]

    si_c = sample_info[cancer_type == c_type, ]
    d_a  = get_dat_one_cancer(c_type, TCR_a, mut_data, si_c)
    d_b  = get_dat_one_cancer(c_type, TCR_b, mut_data, si_c)
    mut1 = d_a$mut_mat
    if(!(mj %in% rownames(mut1))) next
    ta = d_a$TCR_mat; tb = d_b$TCR_mat
    rownames(ta) = paste0("alpha_", rownames(ta))
    rownames(tb) = paste0("beta_",  rownames(tb))
    tcr1 = rbind(ta, tb)

    sel_a = r_a[mut == mj & pval_LRT < pcA, TCR]
    sel_b = r_b[mut == mj & pval_LRT < pcB, TCR]
    idx   = which(rownames(tcr1) %in% c(sel_a, sel_b))
    if(length(idx) == 0) next
    pred = colSums(tcr1[idx, , drop = FALSE] > 0)
    y    = mut1[mj, ]

    sp = si_c[match(colnames(mut1), participant), ]
    acc = rbind(acc, data.table(label = as.integer(y), pred = as.numeric(pred),
          stratum = sp$hot_cold, depth = colSums(tcr1),
          patient = colnames(mut1), cell = paste(c_type, mj)))
  }
  acc[, chain := "both"]
  acc
}

summarise <- function(acc, chain){
  hot = acc[stratum == "hot"]; cold = acc[stratum == "cold"]
  a_hot  = safe_auc(hot$label,  hot$pred)
  a_cold = safe_auc(cold$label, cold$pred)
  ci_hot  = tryCatch(signif(as.numeric(ci.auc(hot$label,  hot$pred,  direction="<", quiet=TRUE)),3), error=function(e) c(NA,NA,NA))
  ci_cold = tryCatch(signif(as.numeric(ci.auc(cold$label, cold$pred, direction="<", quiet=TRUE)),3), error=function(e) c(NA,NA,NA))

  cat(sprintf("\n===== chain = %s : pooled over %d prioritized cells =====\n",
              chain, length(unique(acc$cell))))
  cat(sprintf("  HOT : n=%d (%d mut+),  AUC=%.3f  [%.3f, %.3f]\n",
              nrow(hot), sum(hot$label), a_hot, ci_hot[1], ci_hot[3]))
  cat(sprintf("  COLD: n=%d (%d mut+),  AUC=%.3f  [%.3f, %.3f]\n",
              nrow(cold), sum(cold$label), a_cold, ci_cold[1], ci_cold[3]))
  cat(sprintf("  median TCR depth (recurrent clones): hot=%.0f  cold=%.0f\n",
              median(hot$depth), median(cold$depth)))

  # DeLong test for two AUCs (independent hot vs cold samples)
  rh = roc(hot$label,  hot$pred,  direction="<", quiet=TRUE)
  rc = roc(cold$label, cold$pred, direction="<", quiet=TRUE)
  pv = tryCatch(roc.test(rh, rc, method="delong")$p.value, error=function(e) NA)
  cat(sprintf("  DeLong hot vs cold: p = %.4g\n", pv))

  # depth-matched: rarefy hot patients' detected clones to cold median depth
  target = median(cold$depth)
  # we only have pooled detected-count 'depth' and per-cell predictor; rarefy
  # at the predictor level is not exact, so we approximate the depth effect by
  # restricting hot patients to those whose depth <= hot/cold overlap region:
  ov_hi = quantile(cold$depth, 0.9)
  hot_m = hot[depth <= ov_hi]
  if(length(unique(hot_m$label)) == 2){
    a_hm = safe_auc(hot_m$label, hot_m$pred)
    cat(sprintf("  depth-restricted HOT (depth<=%.0f, cold's 90th pct): n=%d, AUC=%.3f\n",
                ov_hi, nrow(hot_m), a_hm))
  }
  data.table(chain=chain, auc_hot=a_hot, auc_cold=a_cold, n_hot=nrow(hot),
             n_cold=nrow(cold), depth_hot=median(hot$depth),
             depth_cold=median(cold$depth), delong_p=pv)
}

res = list()
for(ch in c("beta","alpha")){
  acc = pool_chain(ch)
  fwrite(acc, file.path(res_dir, sprintf("s6b_pooled_%s.tsv", ch)), sep="\t")
  res[[ch]] = summarise(acc, ch)
}
# combined alpha+beta predictor ("both")
acc_both = pool_both()
fwrite(acc_both, file.path(res_dir, "s6b_pooled_both.tsv"), sep="\t")
res[["both"]] = summarise(acc_both, "both")

cat("\n--- table ---\n"); print(rbindlist(res))

# ---------------------------------------------------------------------
# Write summary tables into the repo (Reviewer comment 3.6, prioritized
# mutations). Reads the pooled patient-level tables just written so the
# block is self-contained / re-runnable.
# ---------------------------------------------------------------------
pooled_summary <- function(ch){
  acc  = fread(file.path(res_dir, sprintf("s6b_pooled_%s.tsv", ch)))
  hot  = acc[stratum == "hot"];  cold = acc[stratum == "cold"]
  cih  = signif(as.numeric(ci.auc(hot$label,  hot$pred,  direction = "<", quiet = TRUE)), 3)
  cic  = signif(as.numeric(ci.auc(cold$label, cold$pred, direction = "<", quiet = TRUE)), 3)
  rh   = roc(hot$label,  hot$pred,  direction = "<", quiet = TRUE)
  rc   = roc(cold$label, cold$pred, direction = "<", quiet = TRUE)
  ovhi = quantile(cold$depth, 0.9); hm = hot[depth <= ovhi]
  data.table(
    chain              = ch,
    prioritized_combos = length(unique(acc$cell)),
    auc_hot            = round(safe_auc(hot$label,  hot$pred),  3),
    hot_95ci           = sprintf("[%.3f, %.3f]", cih[1], cih[3]),
    auc_cold           = round(safe_auc(cold$label, cold$pred), 3),
    cold_95ci          = sprintf("[%.3f, %.3f]", cic[1], cic[3]),
    delong_p           = signif(roc.test(rh, rc, method = "delong")$p.value, 3),
    n_mutpos_hot       = sum(hot$label),
    n_mutpos_cold      = sum(cold$label),
    depth_adj_hot_auc  = round(safe_auc(hm$label, hm$pred), 3))
}
tab_prio = rbindlist(lapply(c("alpha", "beta", "both"), pooled_summary))
fwrite(tab_prio, "s6b_hot_vs_cold_prioritized_summary.txt", sep = "\t")

# prioritized mutation list (Fig 2B/2C set) that was evaluated, incl. combined predictor
prio_list = rbindlist(lapply(c("alpha", "beta", "both"), function(ch){
  combos = fread(file.path(res_dir, sprintf("cv_optimal_pcut_%s_mut_cancer_combos.tsv", ch)))
  data.table(chain = ch, mutation = combos$key_pro, cancer = combos$cancer)
}))
fwrite(prio_list, "s6b_prioritized_mutations.txt", sep = "\t")
cat("wrote: s6b_hot_vs_cold_prioritized_summary.txt, s6b_prioritized_mutations.txt\n")

# ---------------------------------------------------------------------
# Plot: AUC of the prioritized mutations vs. the theoretical AUC null.
# AUC = Mann-Whitney-Wilcoxon statistic; under H0 it is centred at 0.5
# with variance (n1+n0+1)/(12 n1 n0) (Hanley & McNeil, Radiology 1982;
# Mann & Whitney, Ann Math Stat 1947). The representative null SD is
# computed from the prioritized mutations' own carrier counts (n+, n-).
# ---------------------------------------------------------------------
pooled = rbindlist(lapply(c("alpha","beta","both"), function(ch){
  d = fread(file.path(res_dir, sprintf("s6b_pooled_%s.tsv", ch))); d[, chain := ch]; d}))
percell = pooled[, .(np = sum(label == 1), nn = sum(label == 0)), by = .(chain, cell)]
np = median(percell$np); nn = median(percell$nn)
sd0 = sqrt((np + nn + 1) / (12 * np * nn))

prio = rbindlist(lapply(c("alpha","beta","both"), function(ch){
  combos = fread(file.path(res_dir, sprintf("cv_optimal_pcut_%s_mut_cancer_combos.tsv", ch)))
  data.table(chain = ch, auc = combos$auc)}))
h = hist(prio$auc, breaks = seq(0.5, 1, by = 0.025), plot = FALSE)
xx = seq(0, 1, 0.001); nully = dnorm(xx, 0.5, sd0)
ymax2 = max(h$density, nully)
pdf("s6b_auc_hist_prioritized.pdf", width = 6.5, height = 4.5)
plot(h, freq = FALSE, col = "#7570b399", border = "white", xlim = c(0.2, 1),
     ylim = c(0, ymax2*1.05), xlab = "AUC", ylab = "Density",
     main = "Prioritized mutations vs. theoretical AUC null")
lines(xx, nully, col = "black", lwd = 2)
abline(v = 0.5, lty = 2, col = "grey50")
legend("topright", bty = "n", lwd = c(8,2), col = c("#7570b3","black"),
       legend = c(sprintf("prioritized (%d combos, mean %.2f)", nrow(prio), mean(prio$auc)),
                  sprintf("theoretical null N(0.5, %.2f)", sd0)))
dev.off()
cat("wrote: s6b_auc_hist_prioritized.pdf\n")
