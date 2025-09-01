#!/bin/bash

cd /fh/fast/sun_w/yub/tcr/tracerx/tcr-seq

conda activate tcr

# 0. Download file accession list (SraAccList.csv) from SRA
# An example: https://www.ncbi.nlm.nih.gov/sra/SRX5934243[accn]
# alpha chain TCR transcripts running from 5' of reverse transcription reaction to beginning of constant region produced by merging R1 and R2 reads
# Instrument: NextSeq 500
# Strategy: AMPLICON
# Source: TRANSCRIPTOMIC
# Selection: RACE
# Layout: SINGLE


# 1. Fire off all downloads (one Slurm array task per accession)
sbatch -J prefetch -n1 -c1 --mem=1G -t2-0 \
       --array=1-$(wc -l < SraAccList.csv) \
       --wrap='ACC=$(sed -n ${SLURM_ARRAY_TASK_ID}p SraAccList.csv | cut -d, -f1);
               /home/byu2/tools/sratoolkit/sratoolkit.3.2.1-ubuntu64/bin/prefetch \
                 --output-directory raw "$ACC"'


# 2. When the array finishes, verify everything is present (run once in shell)
comm -23 \
  <(cut -d, -f1 SraAccList.csv | sort) \
  <(ls raw | sed "s/\.sra$//" | sort) \
  | tee missing.txt

[[ ! -s missing.txt ]] && echo "✅ All accessions downloaded." || \
echo "⚠️ See missing.txt for accessions that failed."


# 3. Convert sra files to fastq
sbatch \
  --job-name fastq_dump \
  --cpus-per-task 1 \
  --mem 4G \
  --array=1-$(wc -l < SraAccList.csv) \
  --wrap='ID=$(sed -n ${SLURM_ARRAY_TASK_ID}p SraAccList.csv | cut -d, -f1);
          /home/byu2/tools/sratoolkit/sratoolkit.3.2.1-ubuntu64/bin/fastq-dump \
              --split-3 raw/${ID}/${ID}.sra -O fastq/'

# 4. Delete raw data
# rm -r raw
