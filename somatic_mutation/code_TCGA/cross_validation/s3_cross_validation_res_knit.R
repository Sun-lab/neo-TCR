
Sys.setenv(RSTUDIO_PANDOC="/opt/homebrew/bin")

library(rmarkdown)

rmarkdown::render("s3_cross_validation_res.Rmd", 
                  params = list(chain = "alpha"),
                  output_file = "s3_cross_validation_res_alpha.html")
  
rmarkdown::render("s3_cross_validation_res.Rmd", 
                  params = list(chain = "beta"),
                  output_file = "s3_cross_validation_res_beta.html")

sessionInfo()

q(save = "no")

