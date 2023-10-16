library(tidyverse)
library(colorspace)
library(magrittr)
library(kableExtra)


setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/boot_sim_2params_new3_Simulx_SEIR")
source("~/PhD/COVID_France/Dropbox_iris_covid/departement/Données_SPF/Data/data_functions.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")

dir1 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/boot_point_estimates"

# bootstrap replicates
# j = data simulation
# i = bootstrap replicate

SEIR_boot_list_all3 <- list()
point_est_list <- list()
popparams <- c()

# load bootstrap population parameters and calculate bootstrap CIs in 2 ways
for(j in 1:100){
  params_sim_rep <- list()
  
  # load population parameters from individual bootstrap runs 
  for(i in 1:100){
    rm(popparams)
    try(popparams <- read.table(paste0(getwd(), "/SEIR_Simulx_new3_2params_", j, "boot", i, "/populationParameters.txt"), 
                                header = TRUE, sep = ",") %>%
          bootstrap_cleaning(., boot_rep = i, sim_rep = j))
    
    
    if(exists("popparams")){
      params_sim_rep[[i]] <- popparams
    }
  }
  
  # calculate CIs
  pop_params_boot <- bootstrap_CI_calc(params_sim_rep)
  
  
  # load point estimates
  point_est <- read.table(paste0(dir1, "/SEIR_Simulx_new3_pe_", j, "/populationParameters.txt"), 
                                    header = TRUE, sep = ",") %>%
    popparam_cleaning() %>%
    mutate(sim_rep = j) %>%
    rename(mean_est2 = value)
  
  
  SEIR_boot_list_all3[[j]] <- pop_params_boot %>%
    ungroup() %>%
    left_join(point_est, by = c("parameter", "sim_rep")) %>%
    select(parameter, sim_rep, mean_est1, mean_est2, CI_LL1, CI_LL2, CI_UL1, CI_UL2, sd_est)
    
}

save(SEIR_boot_list_all3, file = "SEIR_boot_list_2params_new3_Simulx.RData")



# join point estimates and pivot into longer format 
SEIR_boot_df_long3 <- do.call("rbind.data.frame", SEIR_boot_list_all3) %>%
  pivot_longer(cols = -c(parameter, sim_rep, sd_est), 
               names_to = c(".value", "method"), 
               names_pattern = "(.*)(\\d+)") %>%
  mutate(parameter = factor(parameter, 
                            levels = c("beta_ld1", "beta_BG1"),
                            labels = c("Lockdown", "Barrier gestures")), 
         true_value = ifelse(parameter == "Lockdown", -1.45, -0.8), 
         method = ifelse(method == 1, "Bootstrap SE", "Empirical bootstrap"))


# plot to compare the two CI methods
ggplot(SEIR_boot_df_long3, aes(ymin = CI_LL, ymax = CI_UL, x = sim_rep, y = mean_est, col = method)) + 
  geom_pointrange(position = position_dodge(width = 0.3)) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~parameter, ncol = 1, scale = "free_y") +
  geom_line(aes(y = true_value), linetype = "dashed", col = "darkred", linewidth = 0.8) + 
  labs(title= "Bootstrap SEIR Simulx new3", 
       x = "simulation dataset", y = "coefficient value") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2")


