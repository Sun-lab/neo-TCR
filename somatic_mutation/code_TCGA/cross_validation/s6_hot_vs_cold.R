# =====================================================================
# Does the mutation-prediction performance
# differ between "hot" and "cold" tumors?
#
# Design (mirrors s3_cross_validation_res.Rmd):
#   - Leave-one-cancer-type-out CV; for each (mutation, cancer type) the
#     predictor is the per-patient count of associated TCRs detected, and
#     AUC is computed across patients of that cancer type.
#   - "hot"/"cold" = within-cancer-type median split on T_cell_prop
#     (already in sample_info_intersection.tsv).
#   - For each (mutation, cancer) cell we fix the pcut that maximizes the
#     FULL-cohort AUC (selection is stratum-blind), fix the associated-TCR
#     set, then recompute AUC separately on the hot and cold patients.
#   - Confound check: hot tumors have more TILs -> more TCRs detected ->
#     more signal. We (a) report the TCR-depth gap and (b) re-evaluate hot
#     patients after rarefying their detected-TCR set to the cold depth.
# =====================================================================

suppressMessages({
  library(data.table)
  library(pROC)
  library(reshape2)
})
source("../_utilities.R")

data_dir = "~/research/data/TCR/TCGA/_processed_data"
res_dir  = "~/research/data/TCR/TCGA/_results/_cross_validation"
out_dir  = res_dir

pcuts      = c(5e-3, 1e-2, 5e-2, 1e-1, 2e-1, 5e-1)
MIN_POS    = 2     # min mutation-positive patients required in a stratum
MIN_NEG    = 2     # min mutation-negative patients required in a stratum
N_RAREFY   = 20    # random rarefaction repeats for depth matching
set.seed(2026)

# ---------------------------------------------------------------------
# sample_info + hot/cold assignment (within-cancer median split)
# ---------------------------------------------------------------------
sample_info = fread(file.path(data_dir, "sample_info_intersection.tsv"))
sample_info = sample_info[!is.na(T_cell_prop), ]

sample_info[, med_tcp := median(T_cell_prop), by = cancer_type]
sample_info[, hot_cold := ifelse(T_cell_prop > med_tcp, "hot", "cold")]
cat("hot/cold counts:\n"); print(table(sample_info$hot_cold))

safe_auc <- function(response, predictor){
  if(length(unique(response)) < 2) return(NA_real_)
  as.numeric(auc(response = response, predictor = predictor,
                 quiet = TRUE, direction = "<"))
}

# ---------------------------------------------------------------------
# Core: per (mutation, cancer) AUC_full / AUC_hot / AUC_cold for a chain
# ---------------------------------------------------------------------
eval_chain <- function(chain){
  TCR_data = fread(file.path(data_dir, sprintf("TCR_%s_5.tsv", chain)))
  mut_data = fread(file.path(data_dir, "mutation_freq_geq_5.tsv"))

  res_files = list.files(res_dir, sprintf("train_%s_wo_", chain))
  out = list()

  for(f1 in res_files){
    r1     = fread(file.path(res_dir, f1))
    c_type = sub(sprintf("train_%s_wo_(.*)\\.tsv", chain), "\\1", f1)

    si_c = sample_info[cancer_type == c_type, ]
    if(nrow(si_c) == 0) next

    d1   = get_dat_one_cancer(c_type, TCR_data, mut_data, si_c)
    mut1 = d1$mut_mat
    tcr1 = d1$TCR_mat

    hot_p  = si_c$participant[si_c$hot_cold == "hot"]
    cold_p = si_c$participant[si_c$hot_cold == "cold"]
    is_hot = colnames(mut1) %in% hot_p
    is_cold= colnames(mut1) %in% cold_p

    # per-patient TCR depth = number of distinct detected clones
    depth = colSums(tcr1)

    n_per_mut = rowSums(mut1)
    mut_nms   = rownames(mut1)[n_per_mut >= 3]
    mut_nms   = mut_nms[mut_nms %in% r1$mut]
    if(length(mut_nms) == 0) next

    for(mj in mut_nms){
      y     = mut1[mj, ]
      res_j = r1[mut == mj, ]

      # pick the pcut maximizing full-cohort AUC (stratum-blind selection)
      best = list(auc = -1, pcut = NA, tcr_idx = NULL)
      for(pc in pcuts){
        tcr_sel = res_j$TCR[res_j$pval_LRT < pc]
        idx     = which(rownames(tcr1) %in% tcr_sel)
        if(length(idx) == 0) next
        pred = colSums(tcr1[idx, , drop = FALSE] > 0)
        a    = safe_auc(y, pred)
        if(!is.na(a) && a > best$auc){
          best = list(auc = a, pcut = pc, tcr_idx = idx)
        }
      }
      if(is.null(best$tcr_idx)) next

      pred = colSums(tcr1[best$tcr_idx, , drop = FALSE] > 0)

      yh = y[is_hot];  ph = pred[is_hot]
      yc = y[is_cold]; pc_ = pred[is_cold]

      auc_hot  = if(sum(yh==1) >= MIN_POS && sum(yh==0) >= MIN_NEG) safe_auc(yh, ph) else NA_real_
      auc_cold = if(sum(yc==1) >= MIN_POS && sum(yc==0) >= MIN_NEG) safe_auc(yc, pc_) else NA_real_

      # ---- depth-matched (rarefied) hot AUC --------------------------
      # target depth = median TCR depth among cold patients
      auc_hot_matched = NA_real_
      if(!is.na(auc_hot)){
        target = median(depth[is_cold])
        hot_idx = which(is_hot)
        tcr_rows = best$tcr_idx
        accs = numeric(0)
        for(b in 1:N_RAREFY){
          pred_b = pred
          for(p in hot_idx){
            det = which(tcr1[, p] > 0)            # detected clones in patient p
            if(length(det) > target){
              keep = sample(det, size = target)   # rarefy to cold depth
              pred_b[p] = sum(tcr_rows %in% keep)
            }
          }
          accs = c(accs, safe_auc(yh, pred_b[is_hot]))
        }
        auc_hot_matched = mean(accs, na.rm = TRUE)
      }

      out[[length(out)+1]] = data.table(
        chain = chain, cancer = c_type, mut = mj, pcut = best$pcut,
        n_tcr = length(best$tcr_idx),
        auc_full = best$auc, auc_hot = auc_hot, auc_cold = auc_cold,
        auc_hot_matched = auc_hot_matched,
        n_pos_hot = sum(yh==1), n_pos_cold = sum(yc==1),
        depth_hot = median(depth[is_hot]), depth_cold = median(depth[is_cold]))
    }
    cat(sprintf("  [%s] %s done (%d cells so far)\n", chain, c_type, length(out)))
  }
  rbindlist(out)
}

# ---------------------------------------------------------------------
# Combined alpha+beta predictor ("both"), evaluated over ALL evaluable
# mutations (mirrors eval_chain but selects pcut_A x pcut_B and pools the
# alpha+beta TCR matrices, as in s4_cross_validation_res_both.Rmd).
# ---------------------------------------------------------------------
eval_both <- function(){
  TCR_a    = fread(file.path(data_dir, "TCR_alpha_5.tsv"))
  TCR_b    = fread(file.path(data_dir, "TCR_beta_5.tsv"))
  mut_data = fread(file.path(data_dir, "mutation_freq_geq_5.tsv"))

  res_files = list.files(res_dir, "train_alpha_wo_")
  out = list()

  for(f1 in res_files){
    c_type = sub("train_alpha_wo_(.*)\\.tsv", "\\1", f1)
    r_a = fread(file.path(res_dir, f1))
    r_b = fread(file.path(res_dir, gsub("alpha", "beta", f1)))
    r_a[, TCR := paste0("alpha_", TCR)]
    r_b[, TCR := paste0("beta_",  TCR)]
    r1 = rbind(r_a, r_b)

    si_c = sample_info[cancer_type == c_type, ]
    if(nrow(si_c) == 0) next

    d_a = get_dat_one_cancer(c_type, TCR_a, mut_data, si_c)
    d_b = get_dat_one_cancer(c_type, TCR_b, mut_data, si_c)
    mut1 = d_a$mut_mat
    ta = d_a$TCR_mat; tb = d_b$TCR_mat
    rownames(ta) = paste0("alpha_", rownames(ta))
    rownames(tb) = paste0("beta_",  rownames(tb))
    tcr1 = rbind(ta, tb)

    hot_p  = si_c$participant[si_c$hot_cold == "hot"]
    cold_p = si_c$participant[si_c$hot_cold == "cold"]
    is_hot = colnames(mut1) %in% hot_p
    is_cold= colnames(mut1) %in% cold_p
    depth  = colSums(tcr1)

    n_per_mut = rowSums(mut1)
    mut_nms   = rownames(mut1)[n_per_mut >= 3]
    mut_nms   = mut_nms[mut_nms %in% r1$mut]
    if(length(mut_nms) == 0) next

    for(mj in mut_nms){
      y     = mut1[mj, ]
      res_j = r1[mut == mj, ]

      # pick pcut_A x pcut_B maximizing full-cohort AUC (stratum-blind)
      best = list(auc = -1, pcutA = NA, pcutB = NA, tcr_idx = NULL)
      for(pa in pcuts){
        for(pb in pcuts){
          sel = c(res_j$TCR[res_j$pval_LRT < pa & grepl("^alpha_", res_j$TCR)],
                  res_j$TCR[res_j$pval_LRT < pb & grepl("^beta_",  res_j$TCR)])
          idx = which(rownames(tcr1) %in% sel)
          if(length(idx) == 0) next
          pred = colSums(tcr1[idx, , drop = FALSE] > 0)
          a    = safe_auc(y, pred)
          if(!is.na(a) && a > best$auc){
            best = list(auc = a, pcutA = pa, pcutB = pb, tcr_idx = idx)
          }
        }
      }
      if(is.null(best$tcr_idx)) next

      pred = colSums(tcr1[best$tcr_idx, , drop = FALSE] > 0)
      yh = y[is_hot];  ph = pred[is_hot]
      yc = y[is_cold]; pc_ = pred[is_cold]

      auc_hot  = if(sum(yh==1) >= MIN_POS && sum(yh==0) >= MIN_NEG) safe_auc(yh, ph) else NA_real_
      auc_cold = if(sum(yc==1) >= MIN_POS && sum(yc==0) >= MIN_NEG) safe_auc(yc, pc_) else NA_real_

      auc_hot_matched = NA_real_
      if(!is.na(auc_hot)){
        target = median(depth[is_cold])
        hot_idx = which(is_hot)
        tcr_rows = best$tcr_idx
        accs = numeric(0)
        for(b in 1:N_RAREFY){
          pred_b = pred
          for(p in hot_idx){
            det = which(tcr1[, p] > 0)
            if(length(det) > target){
              keep = sample(det, size = target)
              pred_b[p] = sum(tcr_rows %in% keep)
            }
          }
          accs = c(accs, safe_auc(yh, pred_b[is_hot]))
        }
        auc_hot_matched = mean(accs, na.rm = TRUE)
      }

      out[[length(out)+1]] = data.table(
        chain = "both", cancer = c_type, mut = mj,
        pcut_A = best$pcutA, pcut_B = best$pcutB,
        n_tcr = length(best$tcr_idx),
        auc_full = best$auc, auc_hot = auc_hot, auc_cold = auc_cold,
        auc_hot_matched = auc_hot_matched,
        n_pos_hot = sum(yh==1), n_pos_cold = sum(yc==1),
        depth_hot = median(depth[is_hot]), depth_cold = median(depth[is_cold]))
    }
    cat(sprintf("  [both] %s done (%d cells so far)\n", c_type, length(out)))
  }
  rbindlist(out)
}

# ---------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------
all_res = rbindlist(list(eval_chain("beta"), eval_chain("alpha"), eval_both()),
                    fill = TRUE)
fwrite(all_res, file.path(out_dir, "s6_hot_vs_cold_per_cell.tsv"), sep = "\t")

cat("\n================ SUMMARY ================\n")
for(ch in c("beta","alpha","both")){
  d = all_res[chain == ch & !is.na(auc_hot) & !is.na(auc_cold), ]
  cat(sprintf("\n--- chain = %s : %d paired (mutation x cancer) cells ---\n", ch, nrow(d)))
  if(nrow(d) < 3) next
  cat(sprintf("  mean AUC  hot = %.3f   cold = %.3f   (full = %.3f)\n",
              mean(d$auc_hot), mean(d$auc_cold), mean(d$auc_full)))
  cat(sprintf("  median AUC hot = %.3f   cold = %.3f\n",
              median(d$auc_hot), median(d$auc_cold)))
  w = wilcox.test(d$auc_hot, d$auc_cold, paired = TRUE)
  cat(sprintf("  paired Wilcoxon (hot vs cold): p = %.4g\n", w$p.value))
  cat(sprintf("  median TCR depth  hot = %.0f   cold = %.0f clones\n",
              median(d$depth_hot), median(d$depth_cold)))
  dm = d[!is.na(auc_hot_matched), ]
  if(nrow(dm) >= 3){
    cat(sprintf("  depth-matched: mean AUC hot(rarefied)=%.3f vs cold=%.3f  (paired Wilcoxon p=%.4g)\n",
                mean(dm$auc_hot_matched), mean(dm$auc_cold),
                wilcox.test(dm$auc_hot_matched, dm$auc_cold, paired = TRUE)$p.value))
  }
}
cat("\nwrote: ", file.path(out_dir, "s6_hot_vs_cold_per_cell.tsv"), "\n")

# ---------------------------------------------------------------------
# Write a compact summary table into the repo (
# all evaluable mutations). Reads the per-cell table just written so the
# block is self-contained / re-runnable.
# ---------------------------------------------------------------------
pc_tab = fread(file.path(out_dir, "s6_hot_vs_cold_per_cell.tsv"))
summ_all = rbindlist(lapply(c("alpha", "beta", "both"), function(ch){
  x  = pc_tab[chain == ch & !is.na(auc_hot) & !is.na(auc_cold)]
  dm = x[!is.na(auc_hot_matched)]
  data.table(
    chain             = ch,
    evaluable_cells   = nrow(x),
    mean_auc_hot      = round(mean(x$auc_hot),  3),
    mean_auc_cold     = round(mean(x$auc_cold), 3),
    wilcoxon_p        = signif(wilcox.test(x$auc_hot, x$auc_cold, paired = TRUE)$p.value, 3),
    median_depth_hot  = median(x$depth_hot),
    median_depth_cold = median(x$depth_cold),
    depthmatched_hot  = round(mean(dm$auc_hot_matched), 3),
    depthmatched_cold = round(mean(dm$auc_cold),        3),
    depthmatched_p    = signif(wilcox.test(dm$auc_hot_matched, dm$auc_cold, paired = TRUE)$p.value, 3))
}))
fwrite(summ_all, "s6_hot_vs_cold_summary.txt", sep = "\t")
cat("wrote: s6_hot_vs_cold_summary.txt\n")

# ---------------------------------------------------------------------
# Plot: genome-wide AUC distribution across ALL evaluable mutations,
# computed with the PAPER's prioritization criteria:
#   - alpha alone and beta alone: association p-value cutoff 0.1
#   - combined alpha+beta: the three consensus (pcut_A, pcut_B) pairs
#     {(0.01,0.05),(0.05,0.10),(0.1,0.2)}, taking the best (max) AUC per
#     mutation (matches s5_select_p_threshold.Rmd; for the combined
#     predictor FDR was controlled per (pcut_A,pcut_B) pair in
#     s4_cross_validation_res_both.Rmd).
# AUC equals the Mann-Whitney-Wilcoxon statistic; under H0 it is centred
# at 0.5 with variance (n1+n0+1)/(12 n1 n0) (Hanley & McNeil, Radiology
# 1982; Mann & Whitney, Ann Math Stat 1947). 3 panels, log(1+count),
# x truncated to [0.25,0.90] so the AUC>0.8 signal tail stays visible.
# ---------------------------------------------------------------------
PAIRS = list(c(0.01,0.05), c(0.05,0.10), c(0.10,0.20))   # the three consensus pairs

auc_fixed_chain <- function(chain, cut){
  TCR = fread(file.path(data_dir, sprintf("TCR_%s_5.tsv", chain)))
  mut_data = fread(file.path(data_dir, "mutation_freq_geq_5.tsv"))
  out = list()
  for(f1 in list.files(res_dir, sprintf("train_%s_wo_", chain))){
    r1 = fread(file.path(res_dir, f1))
    c_type = sub(sprintf("train_%s_wo_(.*)\\.tsv", chain), "\\1", f1)
    si_c = sample_info[cancer_type == c_type, ]; if(nrow(si_c)==0) next
    d1 = get_dat_one_cancer(c_type, TCR, mut_data, si_c)
    mut1 = d1$mut_mat; tcr1 = d1$TCR_mat
    mut_nms = rownames(mut1)[rowSums(mut1) >= 3]; mut_nms = mut_nms[mut_nms %in% r1$mut]
    for(mj in mut_nms){
      y = mut1[mj, ]; idx = which(rownames(tcr1) %in% r1[mut==mj & pval_LRT < cut, TCR])
      pred = if(length(idx)==0) rep(0, ncol(tcr1)) else colSums(tcr1[idx,,drop=FALSE] > 0)
      out[[length(out)+1]] = data.table(chain=chain, cancer=c_type, mut=mj,
        auc=safe_auc(y, pred), n_pos=sum(y==1), n_neg=sum(y==0))
    }
  }
  rbindlist(out)
}
# combined predictor: per mutation, MAX AUC over the three consensus pairs
auc_both_merged <- function(pairs){
  TCR_a = fread(file.path(data_dir, "TCR_alpha_5.tsv"))
  TCR_b = fread(file.path(data_dir, "TCR_beta_5.tsv"))
  mut_data = fread(file.path(data_dir, "mutation_freq_geq_5.tsv"))
  out = list()
  for(f1 in list.files(res_dir, "train_alpha_wo_")){
    c_type = sub("train_alpha_wo_(.*)\\.tsv", "\\1", f1)
    r_a = fread(file.path(res_dir, f1)); r_b = fread(file.path(res_dir, gsub("alpha","beta",f1)))
    r_a[, TCR := paste0("alpha_", TCR)]; r_b[, TCR := paste0("beta_", TCR)]
    si_c = sample_info[cancer_type == c_type, ]; if(nrow(si_c)==0) next
    d_a = get_dat_one_cancer(c_type, TCR_a, mut_data, si_c)
    d_b = get_dat_one_cancer(c_type, TCR_b, mut_data, si_c)
    mut1 = d_a$mut_mat
    ta = d_a$TCR_mat; tb = d_b$TCR_mat
    rownames(ta) = paste0("alpha_", rownames(ta)); rownames(tb) = paste0("beta_", rownames(tb))
    tcr1 = rbind(ta, tb)
    mut_nms = rownames(mut1)[rowSums(mut1) >= 3]; mut_nms = mut_nms[mut_nms %in% c(r_a$mut, r_b$mut)]
    for(mj in mut_nms){
      y = mut1[mj, ]
      aucs = sapply(pairs, function(pp){
        idx = which(rownames(tcr1) %in% c(r_a[mut==mj & pval_LRT < pp[1], TCR],
                                          r_b[mut==mj & pval_LRT < pp[2], TCR]))
        pred = if(length(idx)==0) rep(0, ncol(tcr1)) else colSums(tcr1[idx,,drop=FALSE] > 0)
        safe_auc(y, pred) })
      bi = which.max(aucs)
      out[[length(out)+1]] = data.table(chain="both", cancer=c_type, mut=mj, auc=max(aucs),
        best_pair=sprintf("%.2g/%.2g", pairs[[bi]][1], pairs[[bi]][2]),
        auc_pair1=aucs[1], auc_pair2=aucs[2], auc_pair3=aucs[3],
        n_pos=sum(y==1), n_neg=sum(y==0))
    }
  }
  rbindlist(out)
}

allf = rbindlist(list(auc_fixed_chain("alpha", 0.1),
                      auc_fixed_chain("beta",  0.1),
                      auc_both_merged(PAIRS)), fill = TRUE)
fwrite(allf, file.path(out_dir, "s6_auc_paper_criteria_per_cell.tsv"), sep = "\t")

chains_p = c("alpha","beta","both"); cols_p = c("#1b9e77","#d95f02","#7570b3")
titles = c("alpha chain (p < 0.1)", "beta chain (p < 0.1)", "combined alpha+beta (best of 3 pairs)")
brk = seq(-0.005, 1.005, by = 0.01)   # AUC=0.5 is a bin midpoint
XL = c(0.25, 0.90)                     # left truncated; right extended to show AUC>0.8 tail
ymax = max(sapply(chains_p, function(ch){ h = hist(allf[chain==ch & !is.na(auc), auc], breaks=brk, plot=FALSE)
  max(log1p(h$counts[h$mids>=XL[1] & h$mids<=XL[2]])) }))
pdf("s6_auc_density_paper_criteria.pdf", width = 6.5, height = 8)
par(mfrow = c(3,1), mar = c(4,4.2,2.2,1))
for(i in seq_along(chains_p)){
  ch = chains_p[i]; a = allf[chain==ch & !is.na(auc), auc]; h = hist(a, breaks=brk, plot=FALSE)
  plot(h$mids, log1p(h$counts), type="h", lwd=3, lend=1, col=cols_p[i],
       xlim=XL, ylim=c(0, ymax*1.08), xlab="AUC (full-cohort, leave-one-cancer-out CV)",
       ylab="log(1 + count)", main=titles[i])
  abline(v=0.5, lty=2, col="grey40")
  legend("topright", bty="n", cex=0.95,
         legend=c(sprintf("below 0.5:  %d", sum(a<0.5)), sprintf("= 0.5:        %d", sum(abs(a-0.5)<1e-9)),
                  sprintf("above 0.5:  %d", sum(a>0.5)), sprintf("  (AUC>0.8:  %d)", sum(a>0.8)),
                  sprintf("mean AUC = %.3f", mean(a))))
}
dev.off()
cat("wrote: s6_auc_density_paper_criteria.pdf, s6_auc_paper_criteria_per_cell.tsv\n")

# ---------------------------------------------------------------------
# Additional plot: combined alpha+beta AUC distribution at EACH of the
# three consensus (pcut_A, pcut_B) pairs separately (not the best-of-3).
# Each single pair is a clean null (mean ~ 0.5); the looser pairs recover
# the high-AUC (>0.8) signal cells that the strictest pair (0.01/0.05)
# misses. Uses the per-pair AUCs already computed in auc_both_merged().
# ---------------------------------------------------------------------
db = allf[chain == "both"]
pair_cols = c("auc_pair1","auc_pair2","auc_pair3")
pair_lab  = sapply(PAIRS, function(pp) sprintf("alpha p<%.2g, beta p<%.2g", pp[1], pp[2]))
cols3 = c("#3b528b","#21908c","#5dc863")
aucs_list = lapply(pair_cols, function(cl) db[[cl]][!is.na(db[[cl]])])
ymax3 = max(sapply(aucs_list, function(a){ h = hist(a, breaks=brk, plot=FALSE)
  max(log1p(h$counts[h$mids>=XL[1] & h$mids<=XL[2]])) }))
pdf("s6_auc_density_both_3pairs.pdf", width = 6.5, height = 8)
par(mfrow = c(3,1), mar = c(4,4.2,2.2,1))
for(i in 1:3){
  a = aucs_list[[i]]; h = hist(a, breaks=brk, plot=FALSE)
  plot(h$mids, log1p(h$counts), type="h", lwd=3, lend=1, col=cols3[i],
       xlim=XL, ylim=c(0, ymax3*1.08), xlab="AUC (combined alpha+beta, leave-one-cancer-out CV)",
       ylab="log(1 + count)", main=sprintf("combined predictor: %s", pair_lab[i]))
  abline(v=0.5, lty=2, col="grey40")
  legend("topright", bty="n", cex=0.95,
         legend=c(sprintf("cells: %d", length(a)), sprintf("below 0.5:  %d", sum(a<0.5)),
                  sprintf("above 0.5:  %d", sum(a>0.5)), sprintf("AUC>0.8:  %d", sum(a>0.8)),
                  sprintf("mean AUC = %.3f", mean(a))))
}
dev.off()
cat("wrote: s6_auc_density_both_3pairs.pdf\n")
