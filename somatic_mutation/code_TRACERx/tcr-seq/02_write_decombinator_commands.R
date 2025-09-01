setwd('/fh/fast/sun_w/yub/tcr/tracerx/tcr-seq')

library(tidyverse)

# Load meta data
# https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=544699 -> Send to Run Selector
meta <- as_tibble(data.table::fread('SraRunTable.csv'))
meta <- tibble(file_id = meta[['Run']], sample_name = meta[['Sample Name']])

# Create decombinator commands
# Not sure which spacer (I8/M13) to use. Run both and compare to determine
cmd_vec <- plyr::alply(meta, 1, function(x) {
  fn  <- sprintf("%s.fastq", x$file_id)
  ch  <- ifelse(grepl("_TRA", x$sample_name), "a", "b")          # TRA/ TRB → a/ b
  base <- paste(
    "decombinator pipeline",
    "-in",  sprintf("../fastq/%s", fn),
    "-br R2", # In the Chain‑lab 5′‑RACE workflow (Oakes et al., 2017) the barcode is synthesised at the 5′ end and is always read by Read 2
    "-bl 42",
    "-c",  ch,
    "-dz -dc -s -npf"
  )
  c(M13 = paste(base, "-ol M13"),
    I8  = paste(base, "-ol I8"))
})

# write two batch files
dir.create('decombinator_results_M13')
write.table(vapply(cmd_vec, `[`, 1, FUN.VALUE = character(1)),
            "decombinator_results_M13/Batch_M13.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
dir.create('decombinator_results_I8')
write.table(vapply(cmd_vec, `[`, 2, FUN.VALUE = character(1)),
            "decombinator_results_I8/Batch_I8.txt",  quote = FALSE, row.names = FALSE, col.names = FALSE)
