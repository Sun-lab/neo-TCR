
library(data.table)
data_dir = "~/research/data/TCR/TCGA/_processed_data"

sample_ins = fread(file.path(data_dir, "sample_info_intersection.tsv"))

tbl_cts = table(sample_ins[["cancer_type"]])
sort(tbl_cts)

tbl_cts = tbl_cts[tbl_cts > 100]
tbl_cts = as.data.frame(tbl_cts)
tbl_cts

cat(as.character(tbl_cts$Var1), file="cancer_types.tsv", sep="\n")

sessionInfo()

q(save = "no")
