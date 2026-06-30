# =====================================================================
# Candidate neoantigen list (HLA-I). Starting from the mutation x HLA
# enrichment table (step8d), keep (mutation, HLA) pairs where:
#   (1) the HLA allele presents >=1 mutated neo-peptide of the mutation, AND
#   (2) the normalized HLA-allele frequency in carriers (E = carrier freq -
#       background freq) exceeds 0.05  (allele enriched among carriers).
# Two outputs:
#   step8g_neoantigen_list.tsv       -- mutation, HLA, presentation & E columns
#   step8g_neoantigen_list_TCR.tsv   -- same rows + associated TCRs (alpha/beta,
#                                       from the full-TCGA predictors at
#                                       pA=0.1/pB=0.2) and their association
#                                       p-values (pval_LRT).
# =====================================================================
suppressMessages({library(data.table)})
dd="~/research/data/TCR/TCGA"; rd=file.path(dd,"_results")
step8b="/Users/wsun/research/GitHub/neo-TCR/somatic_mutation/code_TCGA/step8b_NetMHCpan"
E_CUT=0.05

## ---- 1. candidate (mutation, HLA) pairs ----
long=fread(file.path(step8b,"step8d_carrier_freq_long.tsv"))
cand=long[class=="HLA-I" & present==TRUE & E>E_CUT]
mm=unique(fread(file.path(rd,"step8_mut_data_HLA_I.tsv"))[,.(key,key_pro)])
cand[,key:=mm$key[match(key_pro,mm$key_pro)]]
setorder(cand,-E)
out=cand[,.(mutation=key_pro, key, HLA, n_neopeptides=P,
            carrier_freq=round(C,4), norm_HLA_freq_in_carriers=round(E,4))]

## ---- 1b. population frequency of the mutation in the TCGA cohort ----
# Numerator: carriers per mutation from mutation_freq_geq_5.tsv (long: participant,key);
# its gene-prefixed key matches `out$key`. Denominator: the same 10,042-patient set
# this file was derived from (sample_info_mutations.tsv, which carries cancer_type).
pd=file.path(dd,"_processed_data")
freq5=fread(file.path(pd,"mutation_freq_geq_5.tsv"))                 # participant, key
sinfo=fread(file.path(pd,"sample_info_mutations.tsv"))[,.(participant,cancer_type)]
N_total=uniqueN(sinfo$participant)                                   # cohort size (10042)
denom=sinfo[,.(n_tot=.N),by=cancer_type]                            # patients per cancer type

car=freq5[key %in% out$key]
car=merge(car, sinfo, by="participant")                              # attach cancer type
# overall: distinct carriers / cohort size
ov=car[,.(n_car=uniqueN(participant)),by=key]
ov[,pop_freq_TCGA:=n_car]      # total carrier count; cohort N_total=10042 (documented in README)
# per cancer type: "n_car/n_tot/CANCER", sorted by carrier count desc
ct=car[,.(n_car=.N),by=.(key,cancer_type)]
ct=merge(ct, denom, by="cancer_type")
setorder(ct,key,-n_car,cancer_type)
ctc=ct[,.(pop_freq_by_cancer=paste(sprintf("%d/%d/%s",n_car,n_tot,cancer_type),collapse=";")),by=key]

out=merge(out, ov[,.(key,pop_freq_TCGA)], by="key", all.x=TRUE, sort=FALSE)
out=merge(out, ctc, by="key", all.x=TRUE, sort=FALSE)
setorder(out,-norm_HLA_freq_in_carriers)

## ---- 1c. presented peptides (sequence, mutation position in peptide, EL %rank) ----
# Re-read the HLA-I NetMHCpan output (same files/logic as step8c) and keep the
# per-peptide presented rows for our candidate (mutation, HLA) pairs. A peptide is
# "presented" when mutant != wild-type, mutant EL %rank < CUT, and mutant < wild-type.
suppressMessages(library(data.table))
net_dir="~/research/data/TCR/TCGA/_results/netMHC"; CUT=2
get_header=function(fnm){
  header=data.frame(fread(fnm,nrows=2,header=TRUE,strip.white=FALSE))
  w=grep("^(HLA|DR)",names(header)); st=unique(w[-1]-w[-length(w)]); st=if(length(st)==1)st else 6
  H=rep("",ncol(header)); for(w1 in w) H[w1:(w1+st-1)]=paste0(names(header)[w1],"_")
  gsub("-","_",paste0(H,header[1,]))}
read_nm=function(f) data.frame(fread(file.path(net_dir,f),col.names=get_header(file.path(net_dir,f)),skip=2))
to_long=function(nm){
  rk=grep("_EL_Rank$",names(nm),value=TRUE)
  m=melt(as.data.table(nm),id.vars=c("ID","Pos","Peptide"),measure.vars=rk,variable.name="HLA",value.name="rank")
  m[,HLA:=gsub("_EL_Rank$","",HLA)]; m[,ID:=gsub("(_MT|_M|_WT|_W)$","",ID)]; m}
# 1-based position of the first residue where the mutant peptide differs from wild-type
first_diff=function(a,b){ av=strsplit(a,"")[[1]]; bv=strsplit(b,"")[[1]]
  n=min(length(av),length(bv)); d=which(av[seq_len(n)]!=bv[seq_len(n)])
  if(length(d)) d[1] else min(length(av),length(bv))+1L }
nm_mut=merge(read_nm("NetMHCpan_4.1b_mut_HLA_A_B.xls"),read_nm("NetMHCpan_4.1b_mut_HLA_C.xls"),by=c("Pos","Peptide","ID"))
nm_ref=merge(read_nm("NetMHCpan_4.1b_ref_HLA_A_B.xls"),read_nm("NetMHCpan_4.1b_ref_HLA_C.xls"),by=c("Pos","Peptide","ID"))
mt=to_long(nm_mut); setnames(mt,c("Peptide","rank"),c("pep_mt","rank_mt"))
rf=to_long(nm_ref); setnames(rf,c("Peptide","rank"),c("pep_wt","rank_wt"))
d=merge(mt,rf,by=c("ID","Pos","HLA"))
d=d[pep_mt!=pep_wt & rank_mt<CUT & rank_mt<rank_wt]            # presented peptides
d[,HLA:=gsub("\\.",":",sub("^HLA[._]","",HLA))]               # "HLA_A01:01" -> "A01:01"
d[,mut_pos:=mapply(first_diff,pep_mt,pep_wt)]
# map NetMHCpan truncated ID -> key_pro (same reconstruction as step8d)
extract_14=function(z){ xx=min(14-nchar(z[1]),nchar(z[2])); paste(z[1],substr(z[2],1,xx),sep="_") }
get_identity=function(x,app){ parts=strsplit(paste(gsub("*","_",x,fixed=TRUE),app,sep="_"),":"); sapply(parts,extract_14) }
idmap=unique(out[,.(key_pro=mutation,HLA)]); idmap[,ID:=gsub("(_MT|_M)$","",get_identity(key_pro,"MT"))]
pep=merge(d[,.(ID,HLA,peptide=pep_mt,mut_pos,EL_rank=round(rank_mt,4))],
          idmap, by=c("ID","HLA"))                            # keep only candidate pairs
setorder(pep,key_pro,HLA,EL_rank)                             # strongest binder first
# widen: peptide_i / mut_pos_i / EL_rank_i  (i = 1..max presented peptides per pair)
pep[,i:=rowid(key_pro,HLA)]
pw=dcast(pep, key_pro+HLA ~ i, value.var=c("peptide","mut_pos","EL_rank"))
ord=unlist(lapply(sort(unique(pep$i)), function(i) paste0(c("peptide","mut_pos","EL_rank"),"_",i)))
setcolorder(pw, c("key_pro","HLA",intersect(ord,names(pw))))
out=merge(out, pw, by.x=c("mutation","HLA"), by.y=c("key_pro","HLA"), all.x=TRUE, sort=FALSE)
setorder(out,-norm_HLA_freq_in_carriers)

fwrite(out, file.path(step8b,"step8g_neoantigen_list.tsv"), sep="\t")
cat(sprintf("candidate (mutation,HLA) pairs (present & E>%.2f): %d  (%d mutations, %d alleles)\n",
            E_CUT, nrow(out), uniqueN(out$mutation), uniqueN(out$HLA)))

## ---- 2. associated TCRs + p-values per mutation ----
TCR_pred=readRDS(file.path(rd,"TCGA_res_TCR_predictors_pval_pA_1.0e-01_pB_2.0e-01.rds"))
aA=fread(file.path(rd,"TCR_alpha_mut_asso_rd_x_with_ct.tsv"), select=c("TCR","mut","pval_LRT"))
aB=fread(file.path(rd,"TCR_beta_mut_asso_rd_x_with_ct.tsv"),  select=c("TCR","mut","pval_LRT"))
setkey(aA,mut,TCR); setkey(aB,mut,TCR)

# format one chain's TCRs as "SEQ(p=..);SEQ(p=..)" and return min p-value
fmt_chain=function(k, seqs, tab){
  if(is.null(seqs) || length(seqs)==0) return(list(str="", minp=NA_real_))
  p=tab[.(k, seqs), pval_LRT]
  list(str=paste(sprintf("%s(p=%.2g)", seqs, p), collapse=";"), minp=suppressWarnings(min(p,na.rm=TRUE)))
}
keys=unique(out$key)
tcr_info=rbindlist(lapply(keys, function(k){
  pr=TCR_pred[[k]]
  fa=fmt_chain(k, pr$TCRA, aA); fb=fmt_chain(k, pr$TCRB, aB)
  data.table(mut_key=k,
             n_TCRA=length(pr$TCRA), n_TCRB=length(pr$TCRB),
             TCRA=fa$str, TCRB=fb$str,
             min_pval_TCRA=fa$minp, min_pval_TCRB=fb$minp)
}))
setnames(tcr_info,"mut_key","key")
out2=merge(out, tcr_info, by="key", all.x=TRUE, sort=FALSE)
out2[, min_TCR_pval := pmin(min_pval_TCRA, min_pval_TCRB, na.rm=TRUE)]
setcolorder(out2, c("mutation","key","HLA","n_neopeptides","carrier_freq",
                    "norm_HLA_freq_in_carriers","n_TCRA","n_TCRB","TCRA","TCRB",
                    "min_pval_TCRA","min_pval_TCRB","min_TCR_pval"))
setorder(out2,-norm_HLA_freq_in_carriers)
fwrite(out2, file.path(step8b,"step8g_neoantigen_list_TCR.tsv"), sep="\t")
cat(sprintf("wrote step8g_neoantigen_list.tsv (%d rows) and step8g_neoantigen_list_TCR.tsv\n", nrow(out)))
cat(sprintf("mutations with associated TCRs available: %d/%d\n",
            sum(keys %in% names(TCR_pred)[!sapply(TCR_pred,is.null)]), length(keys)))
