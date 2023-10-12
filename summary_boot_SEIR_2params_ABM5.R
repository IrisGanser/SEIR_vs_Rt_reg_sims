library(tidyverse)
library(colorspace)
library(magrittr)


setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/boot_sim_2params_ABM5")
source("~/PhD/COVID_France/Dropbox_iris_covid/departement/Données_SPF/Data/data_functions.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")

dir <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/sim_2params_regs"
dir3 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/boot_point_estimates"


# bootstrap replicates
# j = data simulation
# i = bootstrap replicate

ABM5_hybrid_list <- list()
sim_checklist <- list()
popparams <- c()

for(j in 1:100){
  params_sim_rep <- list()
  i_vec <- data.frame(jj = j, ii = rep(NA, 100))
  for(i in 1:100){
    rm(popparams)
    try(popparams <- read.table(paste0(getwd(), "/SEIR_ABM_2params_hybrid5_", j, "boot", i, "/populationParameters.txt"), 
                                header = TRUE, sep = ",") %>%
          bootstrap_cleaning(., boot_rep = i, sim_rep = j))
    
    
    if(exists("popparams")){
      params_sim_rep[[i]] <- popparams
      i_vec[i, 2] <- i
    }
  }
  
  # calculate CIs
  pop_params_boot <- bootstrap_CI_calc(params_sim_rep) %>%
      mutate(model = "hybrid ABM 5")
  
  # load point estimates
  point_est <- read.table(paste0(dir3, "/ABM_hybrid5_pe_", j, "/populationParameters.txt"), 
                          header = TRUE, sep = ",") %>%
    popparam_cleaning() %>%
    mutate(sim_rep = j) %>%
    rename(mean_est2 = value)
  
  
  ABM5_hybrid_list[[j]] <- pop_params_boot %>%
    ungroup() %>%
    left_join(point_est, by = c("parameter", "sim_rep")) %>%
    select(parameter, sim_rep, model, mean_est1, mean_est2, CI_LL1, CI_LL2, CI_UL1, CI_UL2, sd_est)
  
  sim_checklist[[j]] <- i_vec
}

save(ABM5_hybrid_list, file = "ABM5_hybrid_list.RData")

sim_checklist_df <- do.call("rbind.data.frame", sim_checklist) %>% unique()
full_list <- data.frame(jj = rep(1:100, each = 100), 
                        ii = rep(1:100, 100))

diff <- setdiff(full_list, sim_checklist_df)
setdiff(full_list, sim_checklist_df) %>%
  group_by(jj) %>%
  summarize(n = n()) %>%
  print(n = Inf)



# with censored observations
ABM5_rm_list <- list()
sim_checklist <- list()

for(j in 1:100){
  params_sim_rep <- list()
  i_vec <- data.frame(jj = j, ii = rep(NA, 100))
  for(i in 1:100){
    rm(popparams)
    try(popparams <- read.table(paste0(getwd(), "/SEIR_ABM_2params_rm5_", j, "boot", i, "/populationParameters.txt"), 
                                header = TRUE, sep = ",") %>%
          bootstrap_cleaning(., boot_rep = i, sim_rep = j))
    
    
    if(exists("popparams")){
      params_sim_rep[[i]] <- popparams
      i_vec[i, 2] <- i
    }
  }
  
  # calculate CIs
  pop_params_boot <- bootstrap_CI_calc(params_sim_rep) %>%
      mutate(model = "rm ABM 5")
  
  # load point estimates
  point_est <- read.table(paste0(dir3, "/ABM_rm5_pe_", j, "/populationParameters.txt"), 
                          header = TRUE, sep = ",") %>%
    popparam_cleaning() %>%
    mutate(sim_rep = j) %>%
    rename(mean_est2 = value)
  
  ABM5_rm_list[[j]] <- pop_params_boot %>%
    ungroup() %>%
    left_join(point_est, by = c("parameter", "sim_rep")) %>%
    select(parameter, sim_rep, model, mean_est1, mean_est2, CI_LL1, CI_LL2, CI_UL1, CI_UL2, sd_est)
  
  sim_checklist[[j]] <- i_vec
}

save(ABM5_rm_list, file = "ABM5_rm_list.RData")

sim_checklist_df <- do.call("rbind.data.frame", sim_checklist) %>% unique()
full_list <- data.frame(jj = rep(1:100, each = 100), 
                        ii = rep(1:100, 100))

diff <- setdiff(full_list, sim_checklist_df)
setdiff(full_list, sim_checklist_df) %>%
  group_by(jj) %>%
  summarize(n = n()) %>%
  print(n = Inf)


        # summary
load("ABM5_hybrid_list.RData")
load("ABM5_rm_list.RData")

ABM5_boot_df_comp <- bootstrap_summary(ABM5_rm_list) %>% 
  mutate(model = "ABM rm 5") %>%
  bind_rows(bootstrap_summary(ABM5_hybrid_list) %>% mutate(model = "ABM hybrid 5"))

ggplot(ABM5_boot_df_comp, aes(ymin = CI_LL2, ymax = CI_UL2, x = sim_rep, y = mean_est2, col = model)) + 
  geom_pointrange(position = position_dodge(width = 0.5)) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~parameter, ncol = 1, scale = "free_y") +
  geom_line(aes(y = true_value), linetype = "dashed", col = "darkred", linewidth = 0.8) + 
  labs(title= "Bootstrap ABM 5", 
       x = "simulation dataset", y = "coefficient value") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2")

