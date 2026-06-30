# =====================================================================
# Is the number of mutations an HLA allele presents (Figure R2-3) related
# to that allele's population frequency? For each tested HLA allele:
#   x = population carrier frequency (fraction of TCGA patients carrying it)
#   y = number of somatic mutations whose neoantigen it presents (step8c)
# HLA-I population freq from OptiType; HLA-II from the TCGA HLA benchmark.
# =====================================================================
suppressMessages({library(data.table); library(ggplot2); library(ggrepel); library(patchwork)})
dd="~/research/data/TCR/TCGA"; rd=file.path(dd,"_results")
step8b="/Users/wsun/research/GitHub/neo-TCR/somatic_mutation/code_TCGA/step8b_NetMHCpan"

gene_cols=c("HLA-A"="#66C2A5","HLA-B"="#FC8D62","HLA-C"="#8DA0CB",
            "HLA-DR"="#E78AC3","HLA-DQ"="#A6D854","HLA-DP"="#FFD92F")
hla_gene=function(h) fifelse(grepl("^A",h),"HLA-A",fifelse(grepl("^B",h),"HLA-B",fifelse(grepl("^C",h),"HLA-C",
         fifelse(grepl("^DRB",h),"HLA-DR",fifelse(grepl("DQ",h),"HLA-DQ","HLA-DP")))))

## ---- HLA-I: population carrier frequency (OptiType) + n_mut ----
opt=fread(file.path(dd,"OptiTypeCallsHLA_20171207.tsv")); opt[,participant:=substr(aliquot_id,1,12)]
setorder(opt,-Objective); opt=opt[!duplicated(participant)]
al=melt(opt, id.vars="participant", measure.vars=c("A1","A2","B1","B2","C1","C2"),
        value.name="allele")[,.(participant, allele=gsub("\\*","",allele))]
al=unique(al[allele!="" & !is.na(allele)]); bgI=uniqueN(al$participant)
pf_I=al[,.(pop_freq=uniqueN(participant)/bgI), by=allele]
mI=fread(file.path(step8b,"step8c_mut_by_HLA_I_matrix.tsv")); hkI=setdiff(names(mI),"mutation")
nmI=data.table(raw=hkI, allele=gsub("\\.",":",sub("^HLA\\.","",hkI)),
               n_mut=colSums(as.matrix(mI[,..hkI])>0))
dI=merge(nmI, pf_I, by="allele", all.x=TRUE); dI[,class:="HLA-I"]

## ---- HLA-II: population carrier frequency (benchmark) + n_mut ----
suppressMessages(library(readxl))
bm=as.data.table(read_excel("~/research/GitHub/_others/TCGA_HLA_benchmark/HLA_allele_data/TCGA_HLA_benchmark_20200810.xlsx",
                            sheet="HLA_genotyping"))
setnames(bm,"HLA-ClassII","HLA_ClassII"); bm=bm[!is.na(HLA_ClassII)]
ii=bm[,.(allele=unlist(strsplit(HLA_ClassII,"\\|"))), by=.(participant=Sample_Barcode)]
ii[,norm:=gsub("[*:]","",allele)]; bg2=uniqueN(ii$participant)
pts_comp=split(ii$participant, ii$norm); bm_loci=unique(substr(ii$norm,1,4))
parse_ii=function(col){ if(grepl("^DRB",col)) return(gsub("_","",col)); strsplit(sub("^HLA\\.","",col),"\\.")[[1]] }
clean_ii=function(col){ if(grepl("^DRB",col)) return(sub("_","*",col)); paste(strsplit(sub("^HLA\\.","",col),"\\.")[[1]],collapse="-") }
mII=fread(file.path(step8b,"step8c_mut_by_HLA_II_matrix.tsv")); hkII=setdiff(names(mII),"mutation")
popf_II=sapply(hkII,function(col){ comps=parse_ii(col); comps=comps[substr(comps,1,4)%in%bm_loci]
  if(length(comps)==0) NA_real_ else length(Reduce(intersect,pts_comp[comps]))/bg2 })
dII=data.table(raw=hkII, allele=sapply(hkII,clean_ii),
               n_mut=colSums(as.matrix(mII[,..hkII])>0), pop_freq=as.numeric(popf_II), class="HLA-II")

## ---- combine, correlate, plot ----
D=rbind(dI[,.(class,allele,n_mut,pop_freq)], dII[,.(class,allele,n_mut,pop_freq)])
D[,gene:=hla_gene(allele)]
D[,lab:=fifelse(class=="HLA-I" | n_mut>=3 | (pop_freq>0.4 & !is.na(pop_freq)), allele, "")]   # label HLA-I, high-count, or high-freq
for(cl in c("HLA-I","HLA-II")){
  x=D[class==cl & !is.na(pop_freq)]
  cat(sprintf("%s: Spearman(pop freq, # mutations presented) = %.3f (p=%.3g, n=%d alleles)\n",
      cl, cor(x$pop_freq,x$n_mut,method="spearman"),
      cor.test(x$pop_freq,x$n_mut,method="spearman")$p.value, nrow(x)))
}
fwrite(D, file.path(step8b,"step8f_popfreq_vs_presentation.tsv"), sep="\t")

mk=function(cl){
  x=D[class==cl & !is.na(pop_freq)]
  ct=suppressWarnings(cor.test(x$pop_freq,x$n_mut,method="spearman"))
  rho=ct$estimate; pv=ct$p.value
  ggplot(x, aes(pop_freq, n_mut, color=gene)) +
    geom_point(size=2.4, alpha=0.9) +
    geom_text_repel(aes(label=lab), size=2.3, max.overlaps=Inf, min.segment.length=0,
                    box.padding=0.4, seed=7, segment.size=0.3, show.legend=FALSE) +
    scale_color_manual(values=gene_cols, name="gene") +
    theme_bw(base_size=11) + ggtitle(sprintf("%s  (Spearman r = %.2f, pval = %.2g)", cl, rho, pv)) +
    theme(plot.title=element_text(size=10.5)) +
    labs(x="population frequency", y="# somatic mutations presented")
}
comb = mk("HLA-I") + mk("HLA-II") + plot_annotation(tag_levels="A", tag_prefix="(", tag_suffix=")")
ggsave(file.path(step8b,"step8f_popfreq_vs_presentation.pdf"), comb, width=9, height=4.3)
suppressMessages(knitr::plot_crop(file.path(step8b,"step8f_popfreq_vs_presentation.pdf")))
cat("wrote step8f_popfreq_vs_presentation.pdf and .tsv\n")
