library(tidyverse)
library(colorspace)
library(magrittr)


setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/boot_sim_2params_ABM10")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")

dir <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/sim_2params_regs"



# bootstrap replicates
# j = data simulation
# i = bootstrap replicate

ABM10_hybrid_list <- list()
sim_checklist <- list()
popparams <- c()

for(j in 1:100){
  params_sim_rep <- list()
  i_vec <- data.frame(jj = j, ii = rep(NA, 100))
  
  for(i in 1:100){
    rm(popparams)
    try(popparams <- read.table(paste0(getwd(), "/SEIR_ABM_2params_hybrid10_", j, "boot", i, "/populationParameters.txt"), 
                                header = TRUE, sep = ",") %>%
          bootstrap_cleaning(., boot_rep = i, sim_rep = j))
    
    
    if(exists("popparams")){
      params_sim_rep[[i]] <- popparams
      i_vec[i, 2] <- i
    }
  }
  
  # calculate CIs
  if(!is_empty(params_sim_rep)){
    pop_params_boot <- bootstrap_CI_calc(params_sim_rep) %>%
        mutate(model = "hybrid ABM 10")
  } 
  
  # load point estimates
  try(point_est <- read.table(paste0(getwd(), "/ABM_hybrid10_pe_", j, "/populationParameters.txt"), 
                          header = TRUE, sep = ",") %>%
    popparam_cleaning() %>%
    mutate(sim_rep = j) %>%
    rename(mean_est2 = value))
  
  
  ABM10_hybrid_list[[j]] <- pop_params_boot %>%
    ungroup() %>%
    left_join(point_est, by = c("parameter", "sim_rep")) %>%
    select(parameter, sim_rep, model, mean_est1, mean_est2, CI_LL1, CI_LL2, CI_UL1, CI_UL2, sd_est)
  
  sim_checklist[[j]] <- i_vec
}

save(ABM10_hybrid_list, file = "ABM10_hybrid_list.RData")

sim_checklist_df <- do.call("rbind.data.frame", sim_checklist) %>% unique()
full_list <- data.frame(jj = rep(1:100, each = 100), 
                        ii = rep(1:100, 100))

diff <- setdiff(full_list, sim_checklist_df)
setdiff(full_list, sim_checklist_df) %>%
  group_by(jj) %>%
  summarize(n = n()) %>%
  print(n = Inf)



# with censored observations
ABM10_rm_list <- list()
sim_checklist <- list()

for(j in 1:100){
  params_sim_rep <- list()
  i_vec <- data.frame(jj = j, ii = rep(NA, 100))
  for(i in 1:100){
    rm(popparams)
    try(popparams <- read.table(paste0(getwd(), "/SEIR_ABM_2params_rm10_", j, "boot", i, "/populationParameters.txt"), 
                                header = TRUE, sep = ",") %>%
          bootstrap_cleaning(., boot_rep = i, sim_rep = j))
    
    
    if(exists("popparams")){
      params_sim_rep[[i]] <- popparams
      i_vec[i, 2] <- i
    }
  }
  
  # calculate CIs
  if(!is_empty(params_sim_rep)){
    pop_params_boot <- bootstrap_CI_calc(params_sim_rep) %>%
        mutate(model = "rm ABM 10")
  }
  
  # load point estimates
  point_est <- read.table(paste0(getwd(), "/ABM_rm10_pe_", j, "/populationParameters.txt"), 
                          header = TRUE, sep = ",") %>%
    popparam_cleaning() %>%
    mutate(sim_rep = j) %>%
    rename(mean_est2 = value)
  
  ABM10_rm_list[[j]] <- pop_params_boot %>%
    ungroup() %>%
    left_join(point_est, by = c("parameter", "sim_rep")) %>%
    select(parameter, sim_rep, model, mean_est1, mean_est2, CI_LL1, CI_LL2, CI_UL1, CI_UL2, sd_est)
  
  sim_checklist[[j]] <- i_vec
}

save(ABM10_rm_list, file = "ABM10_rm_list.RData")

sim_checklist_df <- do.call("rbind.data.frame", sim_checklist) %>% unique()
full_list <- data.frame(jj = rep(1:100, each = 100), 
                        ii = rep(1:100, 100))

diff <- setdiff(full_list, sim_checklist_df)
setdiff(full_list, sim_checklist_df) %>%
  group_by(jj) %>%
  summarize(n = n()) %>%
  print(n = Inf)


        # summary
load("ABM10_hybrid_list.RData")
load("ABM10_rm_list.RData")

ABM10_boot_df_comp <- bootstrap_summary(ABM10_rm_list, true_val_NPI1 = -1.45, true_val_NPI2 = -0.8) %>% 
  bind_rows(bootstrap_summary(ABM10_hybrid_list, true_val_NPI1 = -1.45, true_val_NPI2 = -0.8))

ggplot(ABM10_boot_df_comp, aes(ymin = CI_LL2, ymax = CI_UL2, x = sim_rep, y = mean_est2, col = model)) + 
  geom_pointrange(position = position_dodge(width = 0.5)) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~parameter, ncol = 1, scale = "free_y") +
  geom_line(aes(y = true_value), linetype = "dashed", col = "darkred", linewidth = 0.8) + 
  labs(title= "Bootstrap ABM 10", 
       x = "simulation dataset", y = "coefficient value") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2")



# long format to compare CIs
ABM10_boot_df_comp_long <- ABM10_boot_df_comp %>%
  select(1:9) %>%
  pivot_longer(cols = -c(parameter, sim_rep, model), 
               names_to = c(".value", "method"), 
               names_pattern = "(.*)(\\d+)") %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.8), 
         method = ifelse(method == 1, "Bootstrap SE", "Empirical bootstrap"))

ggplot(ABM10_boot_df_comp_long, aes(ymin = CI_LL, ymax = CI_UL, x = sim_rep, y = mean_est, col = method)) + 
  geom_pointrange(position = position_dodge(width = 0.5)) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(rows = vars(parameter), cols = vars(model), scales = "free_y") +
  geom_line(aes(y = true_value), linetype = "dashed", col = "darkred", linewidth = 0.8) + 
  labs(title= "Bootstrap ABM 10 Comparison CI methods", 
       x = "simulation dataset", y = "coefficient value") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2")


#### comparison of point estimates ####
point_est_h_list <- list()
point_est_h_nc_list <- list()
point_est_rm_list <- list()
point_est_rm_nc_list <- list()

for(j in 1:5){
  point_est_h <- read.table(paste0(getwd(), "/ABM_hybrid10_pe_", j, "/populationParameters.txt"), 
                          header = TRUE, sep = ",") %>%
    popparam_cleaning() %>%
    mutate(sim_rep = j, model = "ABM 10 hybrid old code") %>%
    rename(mean_est2 = value)
  
  point_est_h_nc <- read.table(paste0(getwd(), "/ABM_hybrid10_pe_nc_", j, "/populationParameters.txt"), 
                            header = TRUE, sep = ",") %>%
    popparam_cleaning() %>%
    mutate(sim_rep = j, model = "ABM 10 hybrid new code") %>%
    rename(mean_est2 = value)
  
  point_est_rm <- read.table(paste0(getwd(), "/ABM_rm10_pe_", j, "/populationParameters.txt"), 
                            header = TRUE, sep = ",") %>%
    popparam_cleaning() %>%
    mutate(sim_rep = j, model = "ABM 10 rm old code") %>%
    rename(mean_est2 = value)
  
  point_est_rm_nc <- read.table(paste0(getwd(), "/ABM_rm10_pe_nc_", j, "/populationParameters.txt"), 
                            header = TRUE, sep = ",") %>%
    popparam_cleaning() %>%
    mutate(sim_rep = j, model = "ABM 10 rm new code") %>%
    rename(mean_est2 = value)
  
  point_est_h_list[[j]] <- point_est_h
  point_est_h_nc_list[[j]] <- point_est_h_nc
  point_est_rm_list[[j]] <- point_est_rm
  point_est_rm_nc_list[[j]] <- point_est_rm_nc
}

point_est_h_df <- do.call("rbind.data.frame", point_est_h_list)
point_est_h_nc_df <- do.call("rbind.data.frame", point_est_h_nc_list)
point_est_rm_df <- do.call("rbind.data.frame", point_est_rm_list)
point_est_rm_nc_df <- do.call("rbind.data.frame", point_est_rm_nc_list)

comp_df_ABM10_nc <- bind_rows(point_est_h_df, point_est_h_nc_df, point_est_rm_df, point_est_rm_nc_df) %>%
  mutate(code = ifelse(grepl("old", model), "old", "new"), 
         model2 = ifelse(grepl("hybrid", model), "hybrid", "random mixing"), 
         parameter = ifelse(grepl("ld", parameter), "NPI 1", "NPI 2"),
         true_value = ifelse(parameter == "NPI 1", -1.45, -0.5))

ggplot(comp_df_ABM10_nc, aes(x = sim_rep, y = mean_est2, col = model2)) +
  geom_point(aes(shape = code)) +
  geom_line(aes(y = true_value), col = "darkred", linetype = "dashed") + 
  facet_wrap(~parameter) + 
  scale_color_brewer(palette = "Dark2")
