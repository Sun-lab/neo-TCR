# =====================================================================
# Intersect mutation-associated TCRs with cancer-reactive (CR) TCRs and
# export tables for the CR-vs-nonCR TCRdist analysis.
# CR TCRs = unique TRA_cdr3 / TRB_cdr3 from the final cancer-reactive cells
# (CD8_CAT_final.tsv U CD4_CAT_final.tsv). Associated TCRs = full-TCGA
# predictors (pA=0.1/pB=0.2), same set used in the earlier TCRdist analysis.
# =====================================================================
suppressMessages({library(data.table)})
dd="~/research/data/TCR/TCGA"; rd=file.path(dd,"_results")
CATd="~/research/CAT/data"
step8b="/Users/wsun/research/GitHub/neo-TCR/somatic_mutation/code_TCGA/step8b_NetMHCpan"
outdir=file.path(step8b,"step8h_tcrdist")

## CR TCR sets (union of CD8 + CD4 final cancer-reactive cells)
cd8=fread(file.path(CATd,"CD8_CAT_final.tsv")); cd4=fread(file.path(CATd,"CD4_CAT_final.tsv"))
CR_a=unique(c(cd8$TRA_cdr3, cd4$TRA_cdr3)); CR_a=CR_a[CR_a!="" & !is.na(CR_a)]
CR_b=unique(c(cd8$TRB_cdr3, cd4$TRB_cdr3)); CR_b=CR_b[CR_b!="" & !is.na(CR_b)]
cat(sprintf("CR cells: CD8=%d CD4=%d ; unique CR CDR3: alpha=%d beta=%d\n",
            nrow(cd8),nrow(cd4),length(CR_a),length(CR_b)))

## associated TCRs for ALL mutations (pA=0.1/pB=0.2 predictors)
TCR_pred=readRDS(file.path(rd,"TCGA_res_TCR_predictors_pval_pA_1.0e-01_pB_2.0e-01.rds"))
keys=names(TCR_pred)[!sapply(TCR_pred,is.null)]
assoc=rbindlist(lapply(keys,function(k){
  pr=TCR_pred[[k]]
  rbind(
    if(length(pr$TCRA)) data.table(mut=k,chain="alpha",cdr3=pr$TCRA) else NULL,
    if(length(pr$TCRB)) data.table(mut=k,chain="beta", cdr3=pr$TCRB) else NULL)
}))
setnames(assoc,"mut","key")
assoc[,is_CR := (chain=="alpha" & cdr3 %in% CR_a) | (chain=="beta" & cdr3 %in% CR_b)]
fwrite(assoc, file.path(outdir,"step8k_assoc_all.tsv"), sep="\t")

## diagnostics
cat(sprintf("\nmutations with predictors: %d ; associated TCR rows: alpha=%d beta=%d\n",
            length(keys), assoc[chain=="alpha",.N], assoc[chain=="beta",.N]))
cat(sprintf("associated TCRs that are CR: alpha=%d/%d (%.1f%%)  beta=%d/%d (%.1f%%)\n",
            assoc[chain=="alpha" & is_CR,.N], assoc[chain=="alpha",.N], 100*mean(assoc[chain=="alpha"]$is_CR),
            assoc[chain=="beta" & is_CR,.N],  assoc[chain=="beta",.N],  100*mean(assoc[chain=="beta"]$is_CR)))

ov=assoc[,.(n_assoc=.N, n_CR=sum(is_CR), n_nonCR=sum(!is_CR)), by=.(key,chain)]
for(ch in c("alpha","beta")){
  q=ov[chain==ch & n_CR>=5][order(-n_CR)]
  cat(sprintf("\n[%s] mutations with >=5 CR-associated TCRs: %d\n", ch, nrow(q)))
  if(nrow(q)) print(q[,.(key,n_assoc,n_CR,n_nonCR)])
}
## also report within the 25 candidate-neoantigen mutations
cand=unique(fread(file.path(step8b,"step8g_neoantigen_list.tsv"))$key)
cat(sprintf("\nof the %d candidate-neoantigen mutations, with >=5 CR: alpha=%d beta=%d\n",
            length(cand), ov[chain=="alpha" & key%in%cand & n_CR>=5,.N],
            ov[chain=="beta" & key%in%cand & n_CR>=5,.N]))
fwrite(ov, file.path(outdir,"step8k_overlap_counts.tsv"), sep="\t")
cat("\nwrote step8k_assoc_all.tsv, step8k_overlap_counts.tsv\n")
