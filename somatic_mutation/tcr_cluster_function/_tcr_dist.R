# This file is a translation of 
# https://github.com/phbradley/conga/blob/master/conga/tcrdist/tcr_distances.py
# into R version

## the BLOSUM62 scoring matrix
##
## Amino acid substitution matrices from protein blocks.
## Henikoff S, Henikoff JG.
## Proc Natl Acad Sci U S A. 1992 Nov 15;89(22):10915-9.
## PMID: 1438297

## TODO: change the blosum name to an uncommon one

source("_blosum62.R")

gap_character = '.'

AMINO_ACIDS = c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')

bsd4 = list()

for (a in AMINO_ACIDS){
  for (b in AMINO_ACIDS){
    cur_key = paste0(a, "_", b)
    bab = blosum[[cur_key]]
    if (a==b){
      bsd4[[cur_key]] = 0
    }else{ ## different
      stopifnot(bab<4)
      if (bab<0){
        bsd4[[cur_key]] = 4
      }else{
        bsd4[[cur_key]] = 4-bab
      }
    }
  }
}


GAP_PENALTY_V_REGION = 4
GAP_PENALTY_CDR3_REGION = 12 # same as GAP_PENALTY_V_REGION=4 since WEIGHT_CDR3_REGION=3 is not applied
WEIGHT_V_REGION = 1
WEIGHT_CDR3_REGION = 3
DISTANCE_MATRIX = bsd4
ALIGN_CDR3S = FALSE
TRIM_CDR3S = TRUE


sequence_distance_with_gappos <- function(shortseq, longseq, gappos){
  
  if (TRIM_CDR3S){
    ntrim = 3 
    ctrim = 2
  }else{
    ntrim = 0  
    ctrim = 0
  }

  remainder = length(shortseq)-gappos
  dist = 0.0
  
  if (ntrim < gappos){
    for (i in seq(ntrim+1, gappos, by=1)){
      dist = dist + DISTANCE_MATRIX[[paste0(shortseq[i], "_", longseq[i])]]
    }
  }
  
  shortseq_rev = rev(shortseq)
  longseq_rev = rev(longseq)
  
  if (ctrim < remainder){
    for (i in seq(ctrim+1, remainder, by=1)){
      dist = dist + DISTANCE_MATRIX[[paste0(shortseq_rev[i], "_", longseq_rev[i])]]
    }
  }
  
  return(dist)
  
}


TcrDistCalculator <- function(chain1, chain2){
  
  seq1 = unlist(strsplit(chain1, ""))
  seq2 = unlist(strsplit(chain2, ""))
  
  if (length(seq1) <= length(seq2)){
    shortseq = seq1
    longseq = seq2
  }else{
    shortseq = seq2
    longseq = seq1    
  }

  lenshort=length(shortseq)
  lenlong=length(longseq)
  lendiff = lenlong - lenshort
  
  stopifnot(lenshort>1)
  stopifnot(lendiff>=0)
  
  if (TRIM_CDR3S){
    # used to say >3+2 but I don't think that's necessary
    stopifnot(lenshort >= 3+2) ## something to align... NOTE: Minimum length of cdr3 protein carried into clones file is currently set in the read_sanger_data.py script!
  }
  
  if (!ALIGN_CDR3S){
    ## if we are not aligning, use a fixed gap position relative to the start of the CDR3
    ## that reflects the typically longer and more variable-length contributions to
    ## the CDR3 from the J than from the V. For a normal-length
    ## CDR3 this would be after the Cys+5 position (ie, gappos = 6; align 6 rsds on N-terminal side of CDR3).
    ## Use an earlier gappos if lenshort is less than 11.
    ##
    gappos = min(6, 3+(lenshort-5)%/%2) # note that gappos in R is the position before gap start
    best_dist = sequence_distance_with_gappos( shortseq, longseq, gappos )
  }
  
  return(WEIGHT_CDR3_REGION * best_dist + lendiff * GAP_PENALTY_CDR3_REGION)
  
}

# chain1="CAAAADYKLSF"
# chain2="CAAAAGDSNYQLIW"

# TcrDistCalculator(chain1, chain2)