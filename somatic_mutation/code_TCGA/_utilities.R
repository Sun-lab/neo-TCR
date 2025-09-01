
get_TCR_matrix <- function(TCR_data, participants, dichotomize=TRUE){
  
  stopifnot(all(c("Participant", "aaSeq", "abundance") %in% names(TCR_data)))
  TCR_mat = reshape2::acast(TCR_data, aaSeq ~ Participant, 
                            value.var="abundance")
  
  TCR_mat = TCR_mat[,colnames(TCR_mat) %in% participants]
  
  TCR_p_others = setdiff(participants, colnames(TCR_mat))
  extra_dat    = matrix(0, nrow=nrow(TCR_mat), ncol=length(TCR_p_others))
  colnames(extra_dat) = TCR_p_others
  TCR_mat = cbind(TCR_mat, extra_dat)
  
  TCR_mat[is.na(TCR_mat)] = 0
  
  if(dichotomize){ 
    TCR_mat[TCR_mat > 0] = 1
  }
  
  return (TCR_mat)
}


get_mut_matrix <- function(mut_data, participants){
  
  stopifnot(all(c("participant", "key") %in% names(mut_data)))
  mut_data$value = 1
  mut_mat = reshape2::acast(mut_data, key ~ participant, 
                            value.var = "value", fill=0)
  
  mut_mat = mut_mat[,colnames(mut_mat) %in% participants]
  
  mut_p_others = setdiff(participants, colnames(mut_mat))
  extra_dat    = matrix(0, nrow=nrow(mut_mat), ncol=length(mut_p_others))
  colnames(extra_dat) = mut_p_others
  mut_mat  = cbind(mut_mat, extra_dat)
  
  return(mut_mat)
}


get_dat_one_cancer <- function(c_type, TCR_data, mut_data, sample_info){
  # need to use the name of c_type instead of cancer_type, 
  # otherwise data.table think it is cancer_type of sample_info
  test_set = sample_info[sample_info$cancer_type == c_type,]
  stopifnot(nrow(test_set) > 0)
  
#  cat(sprintf("cancer type for testing: %s, n for testing = %d\n", 
#              c_type, nrow(test_set)))
  
  participants = test_set[["participant"]]
  
  TCR_mat      = get_TCR_matrix(TCR_data, participants)
  mut_mat      = get_mut_matrix(mut_data, participants)
  
  dim(TCR_mat)
  dim(mut_mat)
  
  stopifnot(setequal(test_set$participant, colnames(TCR_mat)))
  stopifnot(setequal(test_set$participant, colnames(mut_mat)))
  
  TCR_mat = TCR_mat[,match(test_set$participant, colnames(TCR_mat))]
  mut_mat = mut_mat[,match(test_set$participant, colnames(mut_mat))]
  
  dat = list(TCR_mat = TCR_mat, mut_mat=mut_mat)
}
