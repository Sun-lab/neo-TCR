#!/bin/bash


declare -i i=0

seed_file="ensemble_seeds_20.txt"


while IFS= read -r line
do
  stringarray=($line)
  rseed=${stringarray[0]}
  np_seed=${stringarray[1]}
  tf_seed=${stringarray[2]}

  for feature_choice in both_cdr3 both_cdr3_v alpha_cdr3 beta_cdr3 alpha_cdr3_v beta_cdr3_v; do
    for enc_method in one_hot blosum62 atchley pca; do
      for n_CNN in 1 2; do
        for n_units in [32,16] [64,32] [128,64] [256,128]; do
          for freq_cutoff in 10000; do
            i+=1
            cur_name=${feature_choice}_${enc_method}_${n_CNN}_${n_units}_${freq_cutoff}_${rseed}_${np_seed}_${tf_seed}
            . /home/sliu3/anaconda3/etc/profile.d/conda.sh
            conda activate h_tf
            export OMP_NUM_THREADS=1
            sbatch -t 01-10:00:00 -c 2 --mem=8000 -o sc_TCR_log/sc_TCR_${cur_name}.log \
                  -p campus-new --job-name=${i} \
                  sc_TCR.py --feature_choice ${feature_choice} \
                            --model_dir results/sc_TCR_saved_models_cutoff_${freq_cutoff} \
                            --result_dir results/sc_TCR_results_cutoff_${freq_cutoff} \
                            --score_dir results/sc_TCR_scores_cutoff_${freq_cutoff} \
                            --enc_method ${enc_method} \
                            --freq_cutoff ${freq_cutoff} \
                            --lr 0.0001 \
                            --n_CNN ${n_CNN} \
                            --n_dense 2 \
                            --n_units_str ${n_units} \
                            --dropout_flag True \
                            --p_dropout 0.2 \
                            --rseed ${rseed} \
                            --np_seed ${np_seed} \
                            --tf_seed ${tf_seed}
          done
        done
      done
    done
  done
done < "$seed_file"
