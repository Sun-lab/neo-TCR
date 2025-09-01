setwd('/fh/scratch/delete90/sun_w/yub/tcr/tracerx/rna-seq/')

if (!dir.exists('trust4_results')) dir.create('trust4_results')

library(dplyr)


# Load meta information
meta = as_tibble(data.table::fread('sample_file.csv'))

meta$true_file_name = sapply(strsplit(meta$file_name, '_'), function(x) paste(x[5:length(x)], collapse = '_'))
meta$analysis_id = sapply(strsplit(meta$true_file_name, '\\.'), function(x) x[1])


# Check download
for (i in seq(nrow(meta))) {
  folder = meta$file_accession_id[i]
  file = meta$true_file_name[i]
  
  if (!dir.exists(paste0('raw/', folder))) stop(i)
  if (!file.exists(paste0('raw/', folder, '/', file))) stop(i)
}


# Run TRUST4  
analysis_ids = unique(meta$analysis_id)

commands = parallel::mcmapply(function(id) {
  rows = which(meta$analysis_id == id)
  tmp = sapply(meta[rows, ]$true_file_name, function(x) substr(x, nchar(x) - 11, nchar(x)))
  if (!setequal(tmp, c('fq1.fastq.gz', 'fq2.fastq.gz'))) stop(paste0(id, ': no paired files'))
  
  row1 = rows[tmp == 'fq1.fastq.gz']
  folder1 = meta$file_accession_id[row1]
  file1 = meta$true_file_name[row1]
  fastq1 = paste('raw', folder1, file1, sep = '/')
  
  row2 = rows[tmp == 'fq2.fastq.gz']
  folder2 = meta$file_accession_id[row2]
  file2 = meta$true_file_name[row2]
  fastq2 = paste('raw', folder2, file2, sep = '/')
  
  if (!dir.exists(paste0('trust4_results/', id))) dir.create(paste0('trust4_results/', id))
  dir_out = paste0('trust4_results/', id, '/res')
  
  paste0('TRUST4/run-trust4 -1 ', fastq1, ' -2 ', fastq2, ' -f TRUST4/hg19_bcrtcr.fa --ref TRUST4/human_IMGT+C.fa -o ', dir_out)
}, analysis_ids, SIMPLIFY = T, mc.cores = 10)


# Split commands into 5 batches
cmdl = split(commands, f = cut(seq_along(commands), breaks = 5))

for (i in seq_along(cmdl)) write.table(cmdl[[i]], file = paste0('Batch', i, '.txt'), row.names = F, col.names = F, quote = F)


