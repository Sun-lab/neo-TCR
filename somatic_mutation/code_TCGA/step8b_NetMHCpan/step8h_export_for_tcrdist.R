# =====================================================================
# Export tables for the TCRdist analysis of candidate-neoantigen TCRs.
# For each candidate neoantigen mutation (step8g, HLA-I, E>0.05) export its
# associated TCRs (alpha/beta, full-TCGA predictors pA=0.1/pB=0.2); export
# the eligible control universe of TCRs (all TCRs tested in the association
# screen) with each TCR's carrier count (population-frequency numerator) and
# CDR3 length, used to draw frequency- and length-matched control sets.
# =====================================================================
suppressMessages({library(data.table)})
dd="~/research/data/TCR/TCGA"; rd=file.path(dd,"_results")
step8b="/Users/wsun/research/GitHub/neo-TCR/somatic_mutation/code_TCGA/step8b_NetMHCpan"
outdir=file.path(step8b,"step8h_tcrdist"); dir.create(outdir, showWarnings=FALSE)

## candidate neoantigen mutations (unique keys in step8g list)
lst=fread(file.path(step8b,"step8g_neoantigen_list.tsv"))
keys=unique(lst$key)
mm=unique(fread(file.path(rd,"step8_mut_data_HLA_I.tsv"))[,.(key,key_pro)])

## associated TCRs per mutation/chain
TCR_pred=readRDS(file.path(rd,"TCGA_res_TCR_predictors_pval_pA_1.0e-01_pB_2.0e-01.rds"))
assoc=rbindlist(lapply(keys, function(k){
  pr=TCR_pred[[k]]; if(is.null(pr)) return(NULL)
  kp=mm$key_pro[match(k,mm$key)]
  rbind(
    if(length(pr$TCRA)) data.table(mut=k,key_pro=kp,chain="alpha",cdr3=pr$TCRA) else NULL,
    if(length(pr$TCRB)) data.table(mut=k,key_pro=kp,chain="beta", cdr3=pr$TCRB) else NULL)
}))
setnames(assoc,"mut","key")
fwrite(assoc, file.path(outdir,"assoc_tcrs.tsv"), sep="\t")

## eligible control universe = unique tested TCRs (+ carrier count = TCR_freq)
pool_a=unique(fread(file.path(rd,"TCR_alpha_mut_asso_rd_x_with_ct.tsv"),
                    select=c("TCR","TCR_freq"))[,.(cdr3=TCR,carrier_count=TCR_freq)])
pool_b=unique(fread(file.path(rd,"TCR_beta_mut_asso_rd_x_with_ct.tsv"),
                    select=c("TCR","TCR_freq"))[,.(cdr3=TCR,carrier_count=TCR_freq)])
pool_a[,len:=nchar(cdr3)]; pool_b[,len:=nchar(cdr3)]
fwrite(pool_a, file.path(outdir,"control_pool_alpha.tsv"), sep="\t")
fwrite(pool_b, file.path(outdir,"control_pool_beta.tsv"),  sep="\t")

## diagnostics
ac=assoc[,.N,by=.(key_pro,chain)]
cat(sprintf("candidate neoantigen mutations: %d\n", length(keys)))
cat(sprintf("associated TCRs: alpha=%d beta=%d (unique cdr3: a=%d b=%d)\n",
            assoc[chain=="alpha",.N], assoc[chain=="beta",.N],
            uniqueN(assoc[chain=="alpha"]$cdr3), uniqueN(assoc[chain=="beta"]$cdr3)))
cat("\nset size (>=5 needed) per chain:\n")
print(dcast(ac, key_pro~chain, value.var="N", fill=0)[order(-beta)])
cat(sprintf("\nmutations with >=5 alpha: %d ; >=5 beta: %d\n",
            ac[chain=="alpha" & N>=5, .N], ac[chain=="beta" & N>=5, .N]))
cat(sprintf("\ncontrol pool: alpha=%d beta=%d unique TCRs\n", nrow(pool_a), nrow(pool_b)))
cat("carrier_count quantiles (alpha):\n"); print(quantile(pool_a$carrier_count,c(0,.25,.5,.75,.9,.99,1)))
cat("carrier_count quantiles (beta):\n");  print(quantile(pool_b$carrier_count,c(0,.25,.5,.75,.9,.99,1)))
cat("assoc TCR carrier counts present in pool? alpha:",
    mean(assoc[chain=="alpha"]$cdr3 %in% pool_a$cdr3)*100,"%  beta:",
    mean(assoc[chain=="beta"]$cdr3 %in% pool_b$cdr3)*100,"%\n")
cat("wrote tables to step8h_tcrdist/\n")
