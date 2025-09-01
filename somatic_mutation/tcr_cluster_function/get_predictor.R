
# input: chain, p-value cutoff, height cutoff, somatic mutation, 
#        a vector of CDR3 sequences from each patient
# output: a list of two elements
#        output_list$n_cluster:
#         number of clusters associated with given somatic mutation under
#         given p-value cutoff that are within height cutoff of any CDR3 
#         sequence in given vector of sequences, under given chain type
#        output_list$n_tcr:
#         in given vector of CDR3 sequences, the number of CDR3s 
#         (no requirement on uniqueness) fall within height cutoff of any 
#         of the clusters associated with given somatic mutation under
#         given p-value cutoff, under given chain type

# requirement on input:
#    chain, must be one of "alpha" and "beta"
#    hcut, height cutoff, can only be one of 12 and 20 in integer format
#    pcut, pvalue cutoff, must be in float format
#    somatic mutation must exist in the result matrix of asso under the
#       corresponding pvalue cutoff
#       paths to result matrix of asso follow the format below:
#       sprintf("data/cluster_asso/cluster_asso_%s_%d.tsv", chain, hcut)
#    CDR3 sequences cannot have characters outside of the 20 amino acids
#                   cannot have sequences with length<5

# chain = "alpha"
# hcut = 12
# pcut = 0.01
# mutation = "ACVR1C:2:158395120:G:A"
# cdr3s = c("CAAAADYKLSF", 
#           "CAAAAGDSNYQLIW", 
#           "CAAAAGGTSYGKLTF",
#           "CAAAAGNKLTF", 
#           "CAAAAGTASKLTF", 
#           "CAAAANDYKLSF")

source("./_tcr_dist.R")

library(data.table)

get_predictor_cluster <- function(chain, mutation, cdr3s, 
                                  pcut=0.01, 
                                  hcut=12, 
                                  default_flag=TRUE, # whether to use default cluster asso pvalues
                                  asso_file=NULL, 
                                  member_file=NULL, 
                                  verbose=TRUE){
  
  stopifnot(chain%in%c("alpha", "beta"))
  stopifnot(hcut%in%c(12, 20))
  stopifnot(is.numeric(pcut))
  
  stopifnot(length(cdr3s)>0)
  
  cdr3s_concatenated = paste(cdr3s, collapse = '')
  cdr3s_unique_chars = unique(unlist(strsplit(cdr3s_concatenated, "")))
  
  AMINO_ACIDS = c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                  'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
  stopifnot(all(cdr3s_unique_chars%in%AMINO_ACIDS))
  
  nchars = nchar(cdr3s)
  if (min(nchars)<5){
    stop("Error: some of input CDR3 amino acid sequences have length < 5")
  }
  
  
  if (default_flag){ # load default asso and member files
    
    asso_file = file.path(sprintf("data/cluster_asso/cluster_asso_%s_%d.tsv", 
                                  chain, hcut))
    member_file=file.path(sprintf("data/cluster_member/cluster_member_%s_%d.tsv", 
                                  chain, hcut))
    
  }else{ # if not using default asso files, must specify input files
    
    if ((is.null(asso_file))|(is.null(member_file))){
      stop(paste0("Error: if not using default asso and member files,\n",  
                  "must provide path to asso and member files to use"))
    }
    
    # todo: specify requirement for format inside asso file and member file
  }
  
  df_asso = fread(asso_file)
  df_member = fread(member_file)    
  
  if (!(mutation%in%df_asso$mut)){
    warning(paste0("There is no TCR cluster associated with mutation ", 
                   mutation, " in association results","\nreturn NA as output."))
    return(NA)
  }
  
  
  stopifnot(all(df_asso$cluster%in%df_member$h_cut))
  
  
  res_subset = df_asso[which(df_asso$mut == mutation),,drop = FALSE]
  
  # find clusters with pvalue < pvalue cutoff
  ww = which(res_subset$pval_LRT < pcut)
  
  if(length(ww) == 0){
    if (verbose){
      warning(paste0("There is no TCR cluster associated with mutation ", 
                     mutation, " under", " p-value cutoff ", 
                     as.character(pcut), ",\nreturn 0 as output."))
    }
    return(0)
  }
  
  clusters_sig = res_subset$cluster[ww]
  
  cluster_max_dist_mat = NULL
  
  for (i in 1:length(clusters_sig)){
    
    cur_cluster = clusters_sig[i]
    # find TCRs belonging to this cluster
    df_member_cur = df_member[df_member$h_cut==cur_cluster, ]
    cur_tcrs = df_member_cur$aaseq
    # compute TCR dist matrix between TCRs in this cluster and TCRs from patient
    cur_mat = matrix(NA, nrow=length(cur_tcrs), ncol=length(cdr3s))
    for (j in 1:length(cur_tcrs)){
      chain_1 = cur_tcrs[j]
      for (k in 1:length(cdr3s)){
        chain_2 = cdr3s[k]
        cur_mat[j,k] = TcrDistCalculator(chain_1, chain_2)
      }
    }
    # find the max dist between each TCR from patient and TCRs under this cluster
    cur_mat_max_dist = apply(cur_mat, 2, max)
    cluster_max_dist_mat = rbind(cluster_max_dist_mat, cur_mat_max_dist)
    
  }
  
  # find the min of the distances from each TCR to all clusters
  min_vec = as.numeric(apply(cluster_max_dist_mat, 2, min))
  equal_min_mat = t(t(cluster_max_dist_mat) == min_vec)*1
  min_vec_mask = (min_vec<=hcut)*1
  cluster_match_mat = t(t(equal_min_mat)*min_vec_mask)
  
  output_list = list()

  cur_n_cluster = sum(rowSums(cluster_match_mat)>0)
  cur_n_tcr = sum(min_vec_mask)
  
  output_list[["n_cluster"]] = cur_n_cluster
  output_list[["n_tcr"]] = cur_n_tcr
  
  return(output_list)

}

# output_list = get_predictor_cluster(chain, mutation, cdr3s, pcut, hcut)
# output_list$n_cluster
# output_list$n_tcr
