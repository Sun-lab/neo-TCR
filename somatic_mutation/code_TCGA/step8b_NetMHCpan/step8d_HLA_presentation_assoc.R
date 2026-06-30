# =====================================================================
# Do TCR-associated mutations generate peptides that
# bind the HLA alleles expressed by the patients who CARRY the mutation?
#
# Restricted to the netMHCpan-tested alleles (23 HLA-I, 22 HLA-II).
# For each mutation m and tested allele h:
#   P_m[h] = # mutated neo-peptides of m presented by h (step8c)
#   C_m[h] = fraction of m's carriers (all patients with the mutation)
#            who carry allele h  (OptiType=HLA-I, TCGA benchmark=HLA-II)
#   f0[h]  = frequency of h among all typed patients (background)
#   E_m[h] = C_m[h] - f0[h]   (carrier enrichment, removes pop-frequency)
# Per mutation:
#   effect      = mean(E_m | presenting) - mean(E_m | non-presenting)
#   rho_raw     = Spearman(P_m, C_m)
#   rho_enrich  = Spearman(P_m, E_m)
# Aggregated across mutations (one-sample Wilcoxon + permutation null).
# HLA-I and HLA-II separately. Carriers = ALL patients with the mutation.
# =====================================================================
suppressMessages({library(data.table); library(readxl); library(ggplot2); library(ggrepel)})
set.seed(2026)
data_dir = "~/research/data/TCR/TCGA"
res_dir  = "~/research/data/TCR/TCGA/_results"
step8b   = "/Users/wsun/research/GitHub/neo-TCR/somatic_mutation/code_TCGA/step8b_NetMHCpan"

extract_14 = function(z){ xx = min(14 - nchar(z[1]), nchar(z[2])); paste(z[1], substr(z[2],1,xx), sep="_") }
get_identity <- function(x, appendix){
  parts <- strsplit(paste(gsub("*","_",x,fixed=TRUE), appendix, sep="_"), ":"); sapply(parts, extract_14) }

carr = fread(file.path(data_dir,"_processed_data/mutation_freq_geq_5.tsv"))   # participant, key
carriers_of = split(carr$participant, carr$key)

# ---- per-mutation association given allele->participant sets (restricted to tested alleles) ----
run_assoc <- function(P, hla_keys, mut_map, pts_with, bg_patients, label){
  P[, key_pro := mut_map$key_pro[match(mutation, mut_map$ID)]]
  P[, key     := mut_map$key[match(mutation, mut_map$ID)]]
  f0 = sapply(hla_keys, function(h) length(intersect(pts_with[[h]], bg_patients))/length(bg_patients))
  pm = list(); long = list()
  for(i in 1:nrow(P)){
    k = P$key[i]; if(is.na(k)) next
    carr_i = intersect(carriers_of[[k]], bg_patients); if(length(carr_i) < 3) next
    Pvec = as.numeric(P[i, ..hla_keys]); names(Pvec) = hla_keys
    Cvec = sapply(hla_keys, function(h) length(intersect(pts_with[[h]], carr_i))/length(carr_i))
    Evec = Cvec - f0; pres = Pvec > 0
    if(sum(pres) < 1 || sum(!pres) < 1) next
    pm[[length(pm)+1]] = data.table(key_pro=P$key_pro[i], n_carriers=length(carr_i),
      n_present_HLA=sum(pres),
      eff = mean(Evec[pres]) - mean(Evec[!pres]),
      freq_diff = mean(Cvec[pres]) - mean(Cvec[!pres]),                          # carrier-freq diff (raw)
      rho_raw = suppressWarnings(cor(Pvec,Cvec,method="spearman")),
      rho_enrich = suppressWarnings(cor(Pvec,Evec,method="spearman")),
      rho_bin = suppressWarnings(cor(as.integer(pres),Evec,method="spearman")),  # binary presenting vs normalized carrier freq (E)
      rho_bin_p = suppressWarnings(cor.test(as.integer(pres),Evec,method="spearman")$p.value),
      Evec=list(Evec), pres=list(pres))
    long[[length(long)+1]] = data.table(key_pro=P$key_pro[i], HLA=hla_keys,
      P=Pvec, C=Cvec, E=Evec, present=pres)
  }
  pm = rbindlist(pm); long = rbindlist(long)
  perm_p <- function(d, B=2000){
    if(nrow(d)<3) return(NA); obs = mean(d$eff)
    nv = sapply(1:B, function(b) mean(sapply(1:nrow(d), function(j){
      E=d$Evec[[j]]; pr=sample(d$pres[[j]]); mean(E[pr])-mean(E[!pr]) })))
    (1+sum(nv>=obs))/(B+1)
  }
  cat(sprintf("\n=== %s (n=%d mutations; %d tested alleles) ===\n", label, nrow(pm), length(hla_keys)))
  cat(sprintf("  effect mean=%.4f median=%.4f  %%positive=%.0f%%  Wilcoxon p=%.4g  permutation p=%.4g\n",
      mean(pm$eff), median(pm$eff), 100*mean(pm$eff>0), wilcox.test(pm$eff,mu=0)$p.value, perm_p(pm)))
  cat(sprintf("  Spearman(P, C raw):       mean=%.3f median=%.3f  %%>0=%.0f%%  Wilcoxon p=%.4g\n",
      mean(pm$rho_raw,na.rm=TRUE), median(pm$rho_raw,na.rm=TRUE), 100*mean(pm$rho_raw>0,na.rm=TRUE),
      wilcox.test(pm$rho_raw,mu=0)$p.value))
  cat(sprintf("  Spearman(P, E enrichment): mean=%.3f median=%.3f  %%>0=%.0f%%  Wilcoxon p=%.4g\n",
      mean(pm$rho_enrich,na.rm=TRUE), median(pm$rho_enrich,na.rm=TRUE), 100*mean(pm$rho_enrich>0,na.rm=TRUE),
      wilcox.test(pm$rho_enrich,mu=0)$p.value))
  cat(sprintf("  per-mutation test (presentation vs normalized carrier freq, Spearman): %d/%d p<0.05 (%d positive); BH-FDR<0.05: %d\n",
      sum(pm$rho_bin_p<0.05,na.rm=TRUE), nrow(pm), sum(pm$rho_bin_p<0.05 & pm$rho_bin>0,na.rm=TRUE),
      sum(p.adjust(pm$rho_bin_p,"BH")<0.05,na.rm=TRUE)))
  pm[, class := label]; long[, class := label]
  list(per_mut=pm[, .(class,key_pro,n_carriers,n_present_HLA,eff,freq_diff,rho_raw,rho_enrich,rho_bin,rho_bin_p)], long=long)
}

## ===================== HLA-I (OptiType) =====================
opt = fread(file.path(data_dir,"OptiTypeCallsHLA_20171207.tsv"))
opt[, participant := substr(aliquot_id,1,12)]; setorder(opt,-Objective); opt = opt[!duplicated(participant)]
al = melt(opt, id.vars="participant", measure.vars=c("A1","A2","B1","B2","C1","C2"),
          value.name="allele")[, .(participant, allele=gsub("\\*","",allele))]
al = unique(al[allele!="" & !is.na(allele)])
bg_I = unique(al$participant); pts_with_I = split(al$participant, al$allele)
P_I = fread(file.path(step8b,"step8c_mut_by_HLA_I_matrix.tsv"))
hk_I = setdiff(names(P_I),"mutation"); hk_clean = gsub("\\.",":",sub("^HLA\\.","",hk_I)); setnames(P_I,hk_I,hk_clean)
mm_I = unique(fread(file.path(res_dir,"step8_mut_data_HLA_I.tsv"))[,.(key,key_pro)])
mm_I[, ID := gsub("(_WT|_W)$","",get_identity(key_pro,"WT"))]
resI = run_assoc(P_I, hk_clean, mm_I, pts_with_I, bg_I, "HLA-I")

## ===================== HLA-II (benchmark) =====================
bm = as.data.table(read_excel("~/research/GitHub/_others/TCGA_HLA_benchmark/HLA_allele_data/TCGA_HLA_benchmark_20200810.xlsx",
                              sheet="HLA_genotyping"))
setnames(bm, "HLA-ClassII", "HLA_ClassII"); bm = bm[!is.na(HLA_ClassII)]
ii = bm[, .(allele=unlist(strsplit(HLA_ClassII,"\\|"))), by=.(participant=Sample_Barcode)]
ii[, norm := gsub("[*:]","",allele)]
bg_II = unique(ii$participant); pts_comp = split(ii$participant, ii$norm); bm_loci = unique(substr(ii$norm,1,4))
parse_ii <- function(col){ if(grepl("^DRB",col)) return(gsub("_","",col)); strsplit(sub("^HLA\\.","",col),"\\.")[[1]] }
P_II = fread(file.path(step8b,"step8c_mut_by_HLA_II_matrix.tsv")); hk_II = setdiff(names(P_II),"mutation")
pts_with_II = lapply(hk_II, function(col){ comps=parse_ii(col); comps=comps[substr(comps,1,4)%in%bm_loci]
  if(length(comps)==0) character(0) else Reduce(intersect, pts_comp[comps]) }); names(pts_with_II)=hk_II
mm_II = unique(fread(file.path(res_dir,"step8_mut_data_HLA_II.tsv"))[,.(key,key_pro)])
mm_II[, ID := gsub("(_WT|_W)$","",get_identity(key_pro,"WT"))]
resII = run_assoc(P_II, hk_II, mm_II, pts_with_II, bg_II, "HLA-II")

## ===================== save tables =====================
per_mut = rbind(resI$per_mut, resII$per_mut)
long    = rbind(resI$long,    resII$long)
fwrite(per_mut, file.path(step8b,"step8d_assoc_per_mutation.tsv"), sep="\t")
fwrite(long[, .(class,key_pro,HLA,P,C,E,present)], file.path(step8b,"step8d_carrier_freq_long.tsv"), sep="\t")

## ===================== scatter: normalized HLA-freq diff vs Spearman(presenting, carrier freq) =====================
## x = background-controlled enrichment diff (E); presentation as binary; label standout mutations
sc = per_mut[, .(class, key_pro, x=eff, y=rho_bin)]
sc[, lab := ifelse(y > 0.25 | x > 0.05, key_pro, "")]
ps = ggplot(sc, aes(x, y)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey70") +
  geom_vline(xintercept=0, linetype="dashed", color="grey70") +
  geom_point(alpha=0.8, size=2, color="#08519c") +
  geom_text_repel(aes(label=lab), size=2.3, max.overlaps=Inf, min.segment.length=0,
                  box.padding=0.8, point.padding=0.4, force=8, force_pull=0.5, max.iter=40000, seed=7,
                  segment.size=0.3, segment.color="grey60", color="grey20") +
  facet_wrap(~class) + theme_bw(base_size=11) +
  theme(axis.title=element_text(size=10.5)) +
  labs(x="normalized carrier freq. difference (presentation alleles - non-presentation alleles)",
       y=expression(rho*"(neoantigen-presentation, normalized carrier freq.)"))
ggsave(file.path(step8b,"step8d_scatter_eff_vs_spearman.pdf"), ps, width=7, height=4.5)

## ===================== heatmaps: mutation x HLA colored by enrichment E =====================
draw_E_heatmap <- function(d, title, file, hgt){
  m = dcast(d, key_pro ~ HLA, value.var="E"); M = as.matrix(m[,-1]); rownames(M)=m$key_pro
  ro = if(nrow(M)>2) rownames(M)[hclust(dist(M))$order] else rownames(M)
  co = if(ncol(M)>2) colnames(M)[hclust(dist(t(M)))$order] else colnames(M)
  d[, key_pro := factor(key_pro, levels=ro)]; d[, HLA := factor(HLA, levels=co)]
  lim = max(abs(d$E), na.rm=TRUE)
  p = ggplot(d, aes(HLA, key_pro)) +
    geom_tile(aes(fill=E), color="grey92", linewidth=0.2) +
    geom_point(data=d[present==TRUE], shape=4, size=1.1, stroke=0.5, color="black") +
    scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b", midpoint=0,
                         limits=c(-lim,lim), name="normalized\nHLA freq\nin carriers") +
    theme_minimal(base_size=9) +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7),
          axis.text.y=element_text(size=6), panel.grid=element_blank()) +
    labs(title=title, subtitle="x = allele presents ≥1 neo-peptide of the mutation",
         x="HLA allele", y="somatic mutation")
  ggsave(file, p, width=7, height=hgt, limitsize=FALSE)
}
draw_E_heatmap(resI$long,  "HLA-I",  "step8d_heatmap_E_HLA_I.pdf",  hgt=8)
draw_E_heatmap(resII$long, "HLA-II", "step8d_heatmap_E_HLA_II.pdf", hgt=6)

cat("\nwrote: step8d_assoc_per_mutation.tsv, step8d_carrier_freq_long.tsv,\n")
cat("       step8d_scatter_eff_vs_spearman.pdf, step8d_heatmap_E_HLA_{I,II}.pdf\n")
