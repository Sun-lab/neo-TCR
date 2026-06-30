# =====================================================================
# Does the TCR -> mutation association survive
# conditioning on HLA, and is it concentrated in HLA-matched carriers?
#
# Each mutation analyzed SEPARATELY (no patient reuse across mutations).
# Associated TCRs: full-TCGA-cohort association at pA=0.1, pB=0.2
#   (TCGA_res_TCR_predictors_pval_pA_1.0e-01_pB_2.0e-01.rds, from step4,
#    built on TCR_{alpha,beta}_mut_asso_rd_x_with_ct.tsv -- NOT CV).
#
# Rank/AUC-based (robust to separation; no model fitting). Per mutation m
# over patients with TCR + HLA-I data:
#   y   = carrier of m;  tcr = # associated TCRs detected;  hla = has a
#         presenting HLA-I allele for m.
#   AUC_tcr      = AUC(tcr -> y)            (TCR->mutation association)
#   AUC_hla      = AUC(hla -> y)            (does HLA alone predict the mutation?)
#   AUC_tcr | matched / unmatched           (conditioning on HLA)
#   dAUC = AUC_matched - AUC_unmatched      (concentration in HLA-matched)
# Aggregate dAUC tested with a label-permutation null (shuffle hla).
# =====================================================================
suppressMessages({library(data.table); library(ggplot2); library(readxl); library(patchwork)})
set.seed(2026)
dd="~/research/data/TCR/TCGA"; pd=file.path(dd,"_processed_data"); rd=file.path(dd,"_results")
step8b="/Users/wsun/research/GitHub/neo-TCR/somatic_mutation/code_TCGA/step8b_NetMHCpan"

fast_auc <- function(y,x){ n1=sum(y==1); n0=sum(y==0); if(n1<1||n0<1) return(NA_real_)
  r=rank(x); (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0) }

## cohort: TCR-typed AND HLA-I typed
si=fread(file.path(pd,"sample_info_intersection.tsv"))
opt=fread(file.path(dd,"OptiTypeCallsHLA_20171207.tsv")); opt[,participant:=substr(aliquot_id,1,12)]
setorder(opt,-Objective); opt=opt[!duplicated(participant)]
parts=intersect(si$participant, opt$participant)
al=melt(opt[participant%in%parts], id.vars="participant",
        measure.vars=c("A1","A2","B1","B2","C1","C2"), value.name="allele")[,.(participant,allele=gsub("\\*","",allele))]
al=unique(al[allele!="" & !is.na(allele)]); pts_with=split(al$participant,al$allele)

carr=fread(file.path(pd,"mutation_freq_geq_5.tsv")); carriers_of=split(carr$participant,carr$key)
TCRa=fread(file.path(pd,"TCR_alpha_5.tsv"))[Participant%in%parts]
TCRb=fread(file.path(pd,"TCR_beta_5.tsv"))[Participant%in%parts]
TCR_pred=readRDS(file.path(rd,"TCGA_res_TCR_predictors_pval_pA_1.0e-01_pB_2.0e-01.rds"))

P=fread(file.path(step8b,"step8c_mut_by_HLA_I_matrix.tsv")); hk=setdiff(names(P),"mutation")
hk_clean=gsub("\\.",":",sub("^HLA\\.","",hk))
mm=unique(fread(file.path(rd,"step8_mut_data_HLA_I.tsv"))[,.(key,key_pro)])
extr=function(z){xx=min(14-nchar(z[1]),nchar(z[2]));paste(z[1],substr(z[2],1,xx),sep="_")}
mm[,ID:=gsub("(_WT|_W)$","",sapply(strsplit(paste(gsub("*","_",key_pro,fixed=TRUE),"WT",sep="_"),":"),extr))]

## build per-mutation y, tcr, hla
M=list()
for(i in 1:nrow(P)){
  kp=mm$key_pro[match(P$mutation[i],mm$ID)]; k=mm$key[match(P$mutation[i],mm$ID)]
  if(is.na(kp) || is.null(TCR_pred[[k]])) next
  y=as.integer(parts%in%carriers_of[[k]]); if(sum(y)<3) next
  nA=TCRa[aaSeq%in%TCR_pred[[k]]$TCRA,.(n=uniqueN(aaSeq)),by=Participant]
  nB=TCRb[aaSeq%in%TCR_pred[[k]]$TCRB,.(n=uniqueN(aaSeq)),by=Participant]
  tcr=setNames(rep(0L,length(parts)),parts)
  tcr[nA$Participant]=tcr[nA$Participant]+nA$n; tcr[nB$Participant]=tcr[nB$Participant]+nB$n
  pres=hk_clean[as.numeric(P[i,..hk])>0]
  hla=as.integer(parts %in% unique(unlist(pts_with[pres])))
  M[[kp]]=list(y=y, tcr=as.integer(tcr), hla=hla, n_carr=sum(y), n_pres=length(pres))
}
cat(sprintf("mutations analyzed: %d (with >=1 presenting HLA-I allele: %d)\n",
            length(M), sum(sapply(M,function(z)z$n_pres>0))))

## per-mutation AUC metrics
res=rbindlist(lapply(names(M), function(kp){ z=M[[kp]]
  mtc=z$hla==1; umc=z$hla==0
  nc_m=sum(z$y[mtc]==1); nc_u=sum(z$y[umc]==1); nn_m=sum(z$y[mtc]==0); nn_u=sum(z$y[umc]==0)
  data.table(key_pro=kp, n_carrier=z$n_carr, n_present_HLA=z$n_pres, n_matched=sum(mtc),
    carr_matched=nc_m, carr_unmatched=nc_u, noncarr_matched=nn_m, noncarr_unmatched=nn_u,
    min_carr=min(nc_m,nc_u), min_noncarr=min(nn_m,nn_u),
    AUC_tcr=fast_auc(z$y,z$tcr), AUC_hla=fast_auc(z$y,z$hla),
    AUC_tcr_matched   = if(nc_m>=3 && nn_m>=3) fast_auc(z$y[mtc],z$tcr[mtc]) else NA_real_,
    AUC_tcr_unmatched = if(nc_u>=3 && nn_u>=3) fast_auc(z$y[umc],z$tcr[umc]) else NA_real_)
}))
res[, dAUC := AUC_tcr_matched - AUC_tcr_unmatched]
fwrite(res, file.path(step8b,"step8e_perMutation_AUC.tsv"), sep="\t")

cat("\n== Conditioning: does HLA alone predict the mutation? ==\n")
cat(sprintf("  AUC(hla->mutation): median=%.3f mean=%.3f  (~0.5 means HLA does NOT predict the mutation)\n",
            median(res$AUC_hla,na.rm=TRUE), mean(res$AUC_hla,na.rm=TRUE)))
cat(sprintf("  AUC(tcr->mutation): median=%.3f  (within HLA-matched=%.3f, within HLA-unmatched=%.3f)\n",
            median(res$AUC_tcr,na.rm=TRUE),
            median(res$AUC_tcr_matched,na.rm=TRUE), median(res$AUC_tcr_unmatched,na.rm=TRUE)))

cat("\n== Concentration in HLA-matched carriers (dAUC = matched - unmatched) ==\n")
d=res[!is.na(dAUC)]
cat(sprintf("  mutations with both strata: %d\n", nrow(d)))
cat(sprintf("  mean dAUC=%.4f median=%.4f  %%positive=%.0f%%\n", mean(d$dAUC), median(d$dAUC), 100*mean(d$dAUC>0)))
cat(sprintf("  sign test (AUC HLA-matched > HLA-unmatched) p=%.4g\n",
    binom.test(sum(d$dAUC>0), nrow(d), 0.5, "greater")$p.value))

## (figures combined with HLA-II below)

## ---- HLA-II version: AUC(TCR->mutation) vs AUC(HLA->mutation) and matched/unmatched ----
## (benchmark HLA-II typing; same TCR predictor; presenting alleles from step8c HLA-II matrix)
bm = as.data.table(read_excel("~/research/GitHub/_others/TCGA_HLA_benchmark/HLA_allele_data/TCGA_HLA_benchmark_20200810.xlsx",
                              sheet="HLA_genotyping"))
setnames(bm,"HLA-ClassII","HLA_ClassII"); bm=bm[!is.na(HLA_ClassII)]
iiV=bm[,.(allele=unlist(strsplit(HLA_ClassII,"\\|"))), by=.(participant=Sample_Barcode)]
iiV[,norm:=gsub("[*:]","",allele)]
bg2=unique(iiV$participant); pts_comp=split(iiV$participant,iiV$norm); bm_loci=unique(substr(iiV$norm,1,4))
parse_ii <- function(col){ if(grepl("^DRB",col)) return(gsub("_","",col)); strsplit(sub("^HLA\\.","",col),"\\.")[[1]] }
parts2=intersect(si$participant, bg2)
TCRa2=fread(file.path(pd,"TCR_alpha_5.tsv"))[Participant%in%parts2]
TCRb2=fread(file.path(pd,"TCR_beta_5.tsv"))[Participant%in%parts2]
P2=fread(file.path(step8b,"step8c_mut_by_HLA_II_matrix.tsv")); hk2=setdiff(names(P2),"mutation")
pts_with2=lapply(hk2,function(col){ comps=parse_ii(col); comps=comps[substr(comps,1,4)%in%bm_loci]
  if(length(comps)==0) character(0) else Reduce(intersect,pts_comp[comps]) }); names(pts_with2)=hk2
mm2=unique(fread(file.path(rd,"step8_mut_data_HLA_II.tsv"))[,.(key,key_pro)])
mm2[,ID:=gsub("(_WT|_W)$","",sapply(strsplit(paste(gsub("*","_",key_pro,fixed=TRUE),"WT",sep="_"),":"),extr))]
res2=rbindlist(lapply(1:nrow(P2),function(i){
  kp=mm2$key_pro[match(P2$mutation[i],mm2$ID)]; k=mm2$key[match(P2$mutation[i],mm2$ID)]
  if(is.na(kp)||is.null(TCR_pred[[k]])) return(NULL)
  y=as.integer(parts2%in%carriers_of[[k]]); if(sum(y)<3) return(NULL)
  nA=TCRa2[aaSeq%in%TCR_pred[[k]]$TCRA,.(n=uniqueN(aaSeq)),by=Participant]
  nB=TCRb2[aaSeq%in%TCR_pred[[k]]$TCRB,.(n=uniqueN(aaSeq)),by=Participant]
  tcr=setNames(rep(0L,length(parts2)),parts2)
  tcr[nA$Participant]=tcr[nA$Participant]+nA$n; tcr[nB$Participant]=tcr[nB$Participant]+nB$n
  pres=hk2[as.numeric(P2[i,..hk2])>0]; hla=as.integer(parts2 %in% unique(unlist(pts_with2[pres])))
  tcri=as.integer(tcr); mtc=hla==1; umc=hla==0
  data.table(key_pro=kp, AUC_tcr=fast_auc(y,tcri), AUC_hla=fast_auc(y,hla),
    AUC_tcr_matched   = if(sum(y[mtc]==1)>=3 && sum(y[mtc]==0)>=3) fast_auc(y[mtc],tcri[mtc]) else NA_real_,
    AUC_tcr_unmatched = if(sum(y[umc]==1)>=3 && sum(y[umc]==0)>=3) fast_auc(y[umc],tcri[umc]) else NA_real_)
}))
res2[, dAUC := AUC_tcr_matched - AUC_tcr_unmatched]
cat(sprintf("HLA-II: %d patients with TCR+HLA-II, %d mutations (matched/unmatched both strata: %d)\n",
            length(parts2), nrow(res2), sum(!is.na(res2$dAUC))))

mkg <- function(a) ggplot(a, aes(AUC_hla, AUC_tcr)) +
  geom_hline(yintercept=0.5,linetype="dashed",color="grey70") +
  geom_vline(xintercept=0.5,linetype="dashed",color="grey70") +
  geom_point(size=2,alpha=0.8,color="#08519c") + theme_bw(base_size=11) +
  labs(x="AUC (presentation HLA -> mutation)", y="AUC (TCR -> mutation)")
comb = (mkg(res) + ggtitle("HLA-I")) + (mkg(res2) + ggtitle("HLA-II")) +
  plot_annotation(tag_levels="A", tag_prefix="(", tag_suffix=")")
ggsave(file.path(step8b,"step8e_AUC_tcr_vs_hla.pdf"), comb, width=6.2, height=3.3)

## combined matched-vs-unmatched scatter (HLA-I and HLA-II)
mkm <- function(dd) ggplot(dd, aes(AUC_tcr_unmatched, AUC_tcr_matched)) +
  geom_abline(slope=1,intercept=0,linetype="dashed",color="grey60") +
  geom_point(size=2,alpha=0.85,color="#08519c") +
  coord_fixed(ratio=1, xlim=c(0.5,1), ylim=c(0.5,1)) + theme_bw(base_size=11) +
  theme(axis.title=element_text(size=9.5)) +
  labs(x="AUC (patients w/o presentation HLAs)", y="AUC (patients w/ presentation HLAs)")
combm = (mkm(d)+ggtitle("HLA-I")) + (mkm(res2[!is.na(dAUC)])+ggtitle("HLA-II")) +
  plot_annotation(tag_levels="A", tag_prefix="(", tag_suffix=")")
ggsave(file.path(step8b,"step8e_AUC_matched_vs_unmatched.pdf"), combm, width=6.2, height=3.6)
suppressMessages(knitr::plot_crop(file.path(step8b,"step8e_AUC_matched_vs_unmatched.pdf")))  # trim whitespace
cat("\nwrote step8e_perMutation_AUC.tsv, step8e_AUC_matched_vs_unmatched.pdf, step8e_AUC_tcr_vs_hla.pdf\n")
