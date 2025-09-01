
Sys.setenv(RSTUDIO_PANDOC="/opt/homebrew/bin")


library(rmarkdown)

rmarkdown::render("step6b_Zheng_et_al_evaluate_exp_TCR.Rmd", 
                  params = list(pcut_A = "0.005", pcut_B = "0.01"),
                  output_file = "step6b_Zheng_et_al_pA_5e-3_pB_1e-2.html")

rmarkdown::render("step6b_Zheng_et_al_evaluate_exp_TCR.Rmd", 
                  params = list(pcut_A = "0.01", pcut_B = "0.05"),
                  output_file = "step6b_Zheng_et_al_pA_1e-2_pB_5e-2.html")

rmarkdown::render("step6b_Zheng_et_al_evaluate_exp_TCR.Rmd", 
                  params = list(pcut_A = "0.05", pcut_B = "0.1"),
                  output_file = "step6b_Zheng_et_al_pA_5e-2_pB_1e-1.html")

rmarkdown::render("step6b_Zheng_et_al_evaluate_exp_TCR.Rmd", 
                  params = list(pcut_A = "0.1", pcut_B = "0.2"),
                  output_file = "step6b_Zheng_et_al_pA_1e-1_pB_2e-1.html")
sessionInfo()

q(save = "no")
