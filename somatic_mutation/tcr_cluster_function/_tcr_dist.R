# This file is a translation of 
# https://github.com/phbradley/conga/blob/master/conga/tcrdist/tcr_distances.py
# into R version

## the BLOSUM62 scoring matrix
##
## Amino acid substitution matrices from protein blocks.
## Henikoff S, Henikoff JG.
## Proc Natl Acad Sci U S A. 1992 Nov 15;89(22):10915-9.
## PMID: 1438297


blosum <- list(
  'S_W'=-3,
  'G_G'=6,
  'E_M'=-2,
  'A_N'=-2,
  'A_Y'=-2,
  'W_Q'=-2,
  'V_N'=-3,
  'F_K'=-3,
  'G_E'=-2,
  'E_D'=2,
  'W_P'=-4,
  'I_T'=-1,
  'F_D'=-3,
  'K_V'=-2,
  'C_Y'=-2,
  'G_D'=-1,
  'T_N'=0,
  'W_W'=11,
  'S_S'=4,
  'K_C'=-3,
  'E_F'=-3,
  'N_L'=-3,
  'A_K'=-1,
  'Q_P'=-1,
  'F_G'=-3,
  'D_S'=0,
  'C_V'=-1,
  'V_T'=0,
  'H_P'=-2,
  'P_V'=-2,
  'I_Q'=-3,
  'F_V'=-1,
  'W_T'=-2,
  'H_F'=-1,
  'P_D'=-1,
  'Q_R'=1,
  'D_Q'=0,
  'K_Q'=1,
  'D_F'=-3,
  'V_W'=-3,
  'T_C'=-1,
  'A_F'=-2,
  'T_H'=-2,
  'A_Q'=-1,
  'Q_T'=-1,
  'V_F'=-1,
  'F_C'=-2,
  'C_R'=-3,
  'V_P'=-2,
  'H_T'=-2,
  'E_L'=-3,
  'F_R'=-3,
  'I_G'=-4,
  'C_Q'=-3,
  'Y_V'=-1,
  'T_A'=0,
  'T_V'=0,
  'Q_V'=-2,
  'S_K'=0,
  'K_K'=5,
  'E_N'=0,
  'N_T'=0,
  'A_H'=-2,
  'A_C'=0,
  'V_S'=-2,
  'Q_H'=0,
  'H_S'=-1,
  'Q_Y'=-1,
  'P_N'=-2,
  'I_Y'=-1,
  'P_G'=-2,
  'F_N'=-3,
  'H_N'=1,
  'K_H'=-1,
  'N_W'=-4,
  'S_Y'=-2,
  'W_N'=-4,
  'D_Y'=-3,
  'E_Q'=2,
  'K_Y'=-2,
  'S_G'=0,
  'Y_S'=-2,
  'G_R'=-2,
  'A_L'=-1,
  'A_G'=0,
  'T_K'=-1,
  'T_P'=-1,
  'M_V'=1,
  'Q_L'=-2,
  'E_S'=0,
  'H_W'=-2,
  'I_D'=-3,
  'K_F'=-3,
  'N_A'=-2,
  'T_I'=-1,
  'Q_N'=0,
  'K_W'=-3,
  'S_C'=-1,
  'Y_Y'=7,
  'G_V'=-3,
  'L_V'=1,
  'A_R'=-1,
  'M_R'=-1,
  'Y_L'=-1,
  'D_C'=-3,
  'P_P'=7,
  'D_H'=-1,
  'Q_Q'=5,
  'I_V'=3,
  'P_F'=-4,
  'I_A'=-1,
  'F_F'=6,
  'K_T'=-1,
  'L_T'=-1,
  'S_Q'=0,
  'W_F'=1,
  'D_A'=-2,
  'E_Y'=-2,
  'K_A'=-1,
  'Q_S'=0,
  'A_D'=-2,
  'L_R'=-2,
  'T_S'=1,
  'A_V'=0,
  'M_N'=-2,
  'Q_D'=0,
  'E_P'=-1,
  'V_V'=4,
  'D_N'=1,
  'I_S'=-2,
  'P_M'=-2,
  'H_D'=-1,
  'I_L'=2,
  'K_N'=0,
  'L_P'=-3,
  'Y_I'=-1,
  'N_I'=-3,
  'T_Q'=-1,
  'Q_F'=-3,
  'S_M'=-1,
  'E_R'=0,
  'Q_W'=-2,
  'G_N'=0,
  'L_Y'=-1,
  'L_N'=-3,
  'A_S'=1,
  'D_T'=-1,
  'S_T'=1,
  'P_S'=-1,
  'V_R'=-3,
  'D_K'=-1,
  'P_H'=-2,
  'H_C'=-3,
  'Q_I'=-3,
  'H_H'=8,
  'I_I'=4,
  'L_W'=-2,
  'L_L'=4,
  'D_R'=-2,
  'S_I'=-2,
  'D_I'=-3,
  'E_A'=-1,
  'K_I'=-3,
  'Q_K'=1,
  'T_D'=-1,
  'A_W'=-3,
  'Y_R'=-2,
  'M_F'=0,
  'S_P'=-1,
  'H_Q'=0,
  'Y_N'=-2,
  'I_P'=-3,
  'E_C'=-4,
  'H_G'=-2,
  'P_E'=-1,
  'Q_M'=0,
  'H_L'=-3,
  'L_S'=-2,
  'L_H'=-3,
  'N_Q'=0,
  'T_Y'=-2,
  'K_G'=-2,
  'S_E'=0,
  'Y_E'=-2,
  'W_R'=-3,
  'V_M'=1,
  'N_R'=0,
  'G_F'=-3,
  'F_Y'=3,
  'L_Q'=-2,
  'M_Y'=-1,
  'A_P'=-1,
  'S_N'=1,
  'C_L'=-1,
  'L_F'=0,
  'D_W'=-4,
  'S_L'=-2,
  'P_R'=-2,
  'P_K'=-1,
  'Y_G'=-3,
  'C_K'=-3,
  'H_K'=-1,
  'Q_A'=-1,
  'I_F'=0,
  'K_D'=-1,
  'N_C'=-3,
  'L_D'=-4,
  'Y_K'=-2,
  'S_A'=1,
  'W_V'=-3,
  'E_I'=-3,
  'V_I'=3,
  'Q_C'=-3,
  'T_G'=-2,
  'T_L'=-1,
  'L_M'=2,
  'A_T'=0,
  'C_H'=-3,
  'P_Y'=-3,
  'S_H'=-1,
  'H_Y'=2,
  'E_K'=1,
  'C_G'=-3,
  'I_C'=-1,
  'Q_E'=2,
  'K_R'=2,
  'T_E'=-1,
  'L_K'=-2,
  'M_W'=-1,
  'N_Y'=-2,
  'N_H'=1,
  'V_E'=-2,
  'Q_G'=-2,
  'Y_D'=-3,
  'F_Q'=-3,
  'G_Y'=-3,
  'L_I'=2,
  'M_Q'=0,
  'R_A'=-1,
  'C_D'=-3,
  'S_V'=-2,
  'D_D'=6,
  'S_D'=0,
  'P_C'=-3,
  'C_C'=9,
  'W_K'=-3,
  'I_N'=-3,
  'K_L'=-2,
  'N_K'=0,
  'L_G'=-4,
  'M_S'=-1,
  'R_C'=-3,
  'R_D'=-2,
  'V_A'=0,
  'W_I'=-3,
  'T_T'=5,
  'F_M'=0,
  'L_E'=-3,
  'M_M'=5,
  'R_E'=0,
  'W_H'=-2,
  'S_R'=-1,
  'E_W'=-3,
  'P_Q'=-1,
  'H_A'=-2,
  'Y_A'=-2,
  'E_H'=0,
  'R_F'=-3,
  'I_K'=-3,
  'N_E'=0,
  'T_M'=-1,
  'T_R'=-1,
  'M_T'=-1,
  'G_S'=0,
  'L_C'=-1,
  'R_G'=-2,
  'Y_M'=-1,
  'N_F'=-3,
  'Y_Q'=-1,
  'N_P'=-2,
  'R_H'=0,
  'W_M'=-1,
  'C_N'=-3,
  'V_L'=1,
  'F_I'=0,
  'G_Q'=-2,
  'L_A'=-1,
  'M_I'=1,
  'R_I'=-3,
  'W_L'=-2,
  'D_G'=-1,
  'D_L'=-4,
  'I_R'=-3,
  'C_M'=-1,
  'H_E'=0,
  'Y_W'=2,
  'G_P'=-2,
  'W_C'=-2,
  'M_P'=-2,
  'N_S'=1,
  'G_W'=-2,
  'M_K'=-1,
  'R_K'=2,
  'D_E'=2,
  'K_E'=1,
  'R_L'=-2,
  'A_I'=-1,
  'V_Y'=-1,
  'W_A'=-3,
  'Y_F'=3,
  'T_W'=-2,
  'V_H'=-3,
  'F_E'=-3,
  'M_E'=-2,
  'R_M'=-1,
  'E_T'=-1,
  'H_R'=0,
  'P_I'=-3,
  'F_T'=-2,
  'C_I'=-1,
  'H_I'=-3,
  'G_T'=-2,
  'I_H'=-3,
  'R_N'=0,
  'C_W'=-2,
  'W_G'=-2,
  'N_M'=-2,
  'M_L'=2,
  'G_K'=-2,
  'M_G'=-3,
  'K_S'=0,
  'E_V'=-2,
  'N_N'=6,
  'V_K'=-2,
  'R_P'=-2,
  'A_M'=-1,
  'W_E'=-3,
  'F_W'=1,
  'C_F'=-2,
  'V_D'=-3,
  'F_A'=-2,
  'G_I'=-4,
  'M_A'=-1,
  'R_Q'=1,
  'C_T'=-1,
  'W_D'=-4,
  'H_V'=-3,
  'S_F'=-2,
  'P_T'=-1,
  'F_P'=-4,
  'C_E'=-4,
  'H_M'=-2,
  'I_E'=-3,
  'G_H'=-2,
  'R_R'=5,
  'K_P'=-1,
  'C_S'=-1,
  'D_V'=-3,
  'M_H'=-2,
  'M_C'=-1,
  'R_S'=-1,
  'D_M'=-3,
  'E_E'=5,
  'K_M'=-1,
  'V_G'=-3,
  'R_T'=-1,
  'A_A'=4,
  'V_Q'=-2,
  'W_Y'=2,
  'F_S'=-2,
  'G_M'=-3,
  'C_P'=-3,
  'E_G'=-2,
  'I_W'=-3,
  'P_A'=-1,
  'F_L'=0,
  'C_A'=0,
  'G_L'=-4,
  'R_V'=-3,
  'T_F'=-2,
  'Y_P'=-3,
  'M_D'=-3,
  'G_C'=-3,
  'R_W'=-3,
  'N_D'=1,
  'N_V'=-3,
  'V_C'=-1,
  'A_E'=-1,
  'Y_H'=2,
  'D_P'=-1,
  'G_A'=0,
  'R_Y'=-2,
  'P_W'=-4,
  'Y_C'=-2,
  'P_L'=-3,
  'F_H'=-1,
  'I_M'=1,
  'Y_T'=-2,
  'N_G'=0,
  'W_S'=-3
)
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