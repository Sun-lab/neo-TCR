Call TCR from TRACERx TCR-seq data. Test somatic mutation prediction models from TCGA data.

## 01_download_and_format_data.sh 
Download TRACERx TCR-seq sra files and convert to fastq files.

## 02_write_decombinator_commands.R
Write commands for slurm jobs running Decombinator to call TCRs from fastq files.

## 03_run_decombinator.sh
Submit Decombinator jobs to call TCRs.

## 04_guess_spacer 
Guess which spacer was used for each sample based on output file size.

## 05_get_tcr
Get CDR3 sequences from Decombinator outputs. 

Output: Sample x TCR matrix.

## 06_compare_stats_with_paper
To make sure of the results, reproduce several plots in their paper.

## 07_validate_TCGA_results_only_tumor
Test TCGA models with only tumor samples.

## 08_validate_TCGA_results_including_blood
Test TCGA models with tumor + blood samples.
