
# install.packages("data.table", lib="/fh/fast/sun_w/Rlibs")
# install.packages("foreach", lib="/fh/fast/sun_w/Rlibs")
# install.packages("doParallel", lib="/fh/fast/sun_w/Rlibs")
# install.packages("doMC", lib="/fh/fast/sun_w/Rlibs")
# install.packages("dplyr", lib="/fh/fast/sun_w/Rlibs")
# install.packages("reshape2", lib="/fh/fast/sun_w/Rlibs")

lib_path = "/fh/fast/sun_w/Rlibs"
data_dir = "/fh/fast/sun_w/research/data/TCR/TCGA/_processed_data"
res_dir  = "/fh/fast/sun_w/research/data/TCR/TCGA/_results/_cross_validation"
at_cluster = TRUE

if(!file.exists(lib_path)){
  lib_path = .libPaths()
  data_dir = "~/research/data/TCR/TCGA/_processed_data"
  res_dir  = "~/research/data/TCR/TCGA/_results/_cross_validation"
  at_cluster = FALSE
}

lib_path

library(data.table, lib.loc=lib_path)
library(foreach, lib.loc=lib_path)
library(doParallel, lib.loc=lib_path)
library(dplyr, lib.loc=lib_path)
library(reshape2, lib.loc=lib_path)
source("../_utilities.R")

args=(commandArgs(TRUE))
args

if (length(args) != 2) {
  message("one argument is expected, use 'BRCA' as default.\n")
  cancer_type = "BRCA"
  chain = "alpha"
}else{
  for(i in 1:2){
    eval(parse(text=args[[i]]))
  }
}

cancer_type
chain

# --------------------------------------------------------------------
# Data analysis
# --------------------------------------------------------------------

## Read in data
TCR_data   <- fread(file.path(data_dir, sprintf("TCR_%s_5.tsv", chain)))
mut_data   <- fread(file.path(data_dir, "mutation_freq_geq_5.tsv"))
sample_ins <- fread(file.path(data_dir, "sample_info_intersection.tsv"))

dim(TCR_data)
dim(mut_data)
dim(sample_ins)

head(TCR_data)
head(mut_data)
head(sample_ins)

ctype = cancer_type
train_set <- sample_ins[sample_ins$cancer_type != ctype,]

cat(sprintf("cancer type for testing: %s, n for training = %d\n", 
            cancer_type, nrow(train_set)))
cat(date(), "\n")

participants <- train_set[["participant"]]

TCR_mat      <- get_TCR_matrix(TCR_data, participants)
mut_mat      <- get_mut_matrix(mut_data, participants)

dim(TCR_mat)
dim(mut_mat)

stopifnot(setequal(train_set$participant, colnames(TCR_mat)))
stopifnot(setequal(train_set$participant, colnames(mut_mat)))

TCR_mat <- TCR_mat[,match(train_set$participant, colnames(TCR_mat))]
mut_mat <- mut_mat[,match(train_set$participant, colnames(mut_mat))]

n_per_mut <- rowSums(mut_mat)
n_per_TCR <- rowSums(TCR_mat)
mut_mat <- mut_mat[which(n_per_mut >= 5), ]
TCR_mat <- TCR_mat[which(n_per_TCR >= 5), ]

dim(TCR_mat)
dim(mut_mat)

cores <- Sys.getenv("SLURM_CPUS_ON_NODE")
cores <- as.numeric(cores)
cores

n_TCRs = nrow(TCR_mat)
n_muts = nrow(mut_mat)

# registerDoParallel(cores=cores)

TCR_rd  = train_set$abundance
TCR_rd2 = colSums(TCR_mat)
summary(TCR_rd)
summary(TCR_rd2)

cor(TCR_rd, TCR_rd2, method="spearman")

TCR_rd = log10(TCR_rd)

logistic_reg <- function(y, x, TCR_rd, l0){
  l1  = glm(y ~ x*TCR_rd, family = "binomial")
  lrt = anova(l0, l1, test = "LRT")
  pv1 = as.numeric(lrt$`Pr(>Chi)`[2])
  pv1
}

res_i = NULL

for(i in 1:n_TCRs){
# res_i = foreach(i = 1:nTCRs, .combine=rbind) %dopar%{
  if(i %% 100 == 0){ cat(i, date(), "\n") }
  y  = TCR_mat[i,]
  l0 = glm(y ~ TCR_rd, family = "binomial")
  
  date()
  pvals <- sapply(1:n_muts, function(j) {
    logistic_reg(y, mut_mat[j,], TCR_rd=TCR_rd, l0=l0)
  })
  date()
  
  df1 = data.frame(TCR = rep(rownames(TCR_mat)[i], n_muts), 
                   mut = rownames(mut_mat), 
                   pval_LRT = pvals)
  df1 = df1[df1$pval_LRT <= 0.5, ]
  # df1
  res_i = rbind(res_i, df1)
}


cat(sprintf("%d TCRs, %d mutations, %d association results\n", 
                 nrow(TCR_mat), nrow(mut_mat), nrow(res_i)))

fwrite(res_i, file=file.path(res_dir, sprintf("train_%s_wo_%s.tsv", 
                                              chain, cancer_type)), 
       sep="\t")

gc()
sessionInfo()

q(save = "no")
