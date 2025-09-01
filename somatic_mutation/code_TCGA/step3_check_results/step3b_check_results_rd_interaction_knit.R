
Sys.setenv(RSTUDIO_PANDOC="/opt/homebrew/bin")


library(rmarkdown)

rmarkdown::render("step3b_check_results_rd_interaction.Rmd", 
                  params = list(chain = "alpha"),
                  output_file = "step3b_check_results_rd_interaction_alpha.html")
  
rmarkdown::render("step3b_check_results_rd_interaction.Rmd", 
                  params = list(chain = "beta"),
                  output_file = "step3b_check_results_rd_interaction_beta.html")

sessionInfo()

q(save = "no")
