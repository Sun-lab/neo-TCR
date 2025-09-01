
# Sys.setenv(RSTUDIO_PANDOC="/app/software/RStudio-Server/1.4.1717-foss-2021b-Java-11-R-4.1.2/bin/pandoc")
# Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")
Sys.setenv(RSTUDIO_PANDOC="/opt/homebrew/bin")


library(rmarkdown)

rmarkdown::render("step2b_TCR_mut_asso_rd_interaction.Rmd", 
                  params = list(chain = "alpha"),
                  output_file = "step2b_TCR_alpha_mut_asso_rd_interaction.html")
  
rmarkdown::render("step2b_TCR_mut_asso_rd_interaction.Rmd", 
                  params = list(chain = "beta"),
                  output_file = "step2b_TCR_beta_mut_asso_rd_interaction.html")

sessionInfo()

q(save = "no")
