library(tidyverse)
library(colorspace)
library(magrittr)
library(kableExtra)


setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/boot_sim_2params_new4_Simulx_SEIR")
source("~/PhD/COVID_France/Dropbox_iris_covid/departement/Données_SPF/Data/data_functions.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")



# bootstrap replicates
# j = data simulation
# i = bootstrap replicate

SEIR_boot_list_all4 <- list()
point_est_list <- list()
popparams <- c()

# load bootstrap population parameters and calculate bootstrap CIs in 2 ways
for(j in 1:50){
  params_sim_rep <- list()
  
  # load population parameters from individual bootstrap runs 
  for(i in 1:100){
    rm(popparams)
    try(popparams <- read.table(paste0(getwd(), "/SEIR_Simulx_new4_2params_", j, "boot", i, "/populationParameters.txt"), 
                                header = TRUE, sep = ",") %>%
          bootstrap_cleaning(., boot_rep = i, sim_rep = j))
    
    
    if(exists("popparams")){
      params_sim_rep[[i]] <- popparams
    }
  }
  
  # calculate CIs
  if(!is_empty(params_sim_rep)){
    pop_params_boot <- bootstrap_CI_calc(params_sim_rep)
  }
  
  
  # load point estimates
  point_est <- read.table(paste0(getwd(), "/SEIR_Simulx_new4_pe_", j, "/populationParameters.txt"), 
                                    header = TRUE, sep = ",") %>%
    popparam_cleaning() %>%
    mutate(sim_rep = j) %>%
    rename(mean_est2 = value)
  
  
  SEIR_boot_list_all4[[j]] <- pop_params_boot %>%
    ungroup() %>%
    left_join(point_est, by = c("parameter", "sim_rep")) %>%
    select(parameter, sim_rep, mean_est1, mean_est2, CI_LL1, CI_LL2, CI_UL1, CI_UL2, sd_est)
    
}

save(SEIR_boot_list_all4, file = "SEIR_boot_list_2params_new4_Simulx.RData")



# join point estimates and pivot into longer format 
SEIR_boot_df_long4 <- do.call("rbind.data.frame", SEIR_boot_list_all4) %>%
  pivot_longer(cols = -c(parameter, sim_rep, sd_est), 
               names_to = c(".value", "method"), 
               names_pattern = "(.*)(\\d+)") %>%
  mutate(parameter = factor(parameter, 
                            levels = c("beta_ld1", "beta_BG1"),
                            labels = c("Lockdown", "Barrier gestures")), 
         true_value = ifelse(parameter == "Lockdown", -1.45, -0.5), 
         method = ifelse(method == 1, "Bootstrap SE", "Empirical bootstrap"))


# plot to compare the two CI methods
ggplot(SEIR_boot_df_long4, aes(ymin = CI_LL, ymax = CI_UL, x = sim_rep, y = mean_est, col = method)) + 
  geom_pointrange(position = position_dodge(width = 0.3)) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~parameter, ncol = 1, scale = "free_y") +
  geom_line(aes(y = true_value), linetype = "dashed", col = "darkred", linewidth = 0.8) + 
  labs(title= "Bootstrap SEIR Simulx new4", 
       x = "simulation dataset", y = "coefficient value") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2")



## just point estimates
point_est_list <- list()
point_est_short_list <- list()
point_est_short_initE_list <- list()
point_est_initE_list <- list()
for(j in 1:100){
  point_est <- read.table(paste0(getwd(), "/SEIR_Simulx_new4_pe_", j, "/populationParameters.txt"), 
                          header = TRUE, sep = ",") %>%
    popparam_cleaning() %>%
    mutate(sim_rep = j) %>%
    rename(mean_est2 = value)
  
  point_est_list[[j]] <- point_est
  
  
  point_est_short <- read.table(paste0(getwd(), "/SEIR_Simulx_new4_pe_short_", j, "/populationParameters.txt"), 
                          header = TRUE, sep = ",") %>%
    popparam_cleaning() %>%
    mutate(sim_rep = j) %>%
    rename(mean_est2 = value)
  
  point_est_short_list[[j]] <- point_est_short
  
  
  point_est_initE <- read.table(paste0(getwd(), "/SEIR_Simulx_new4_pe_initE", j, "/populationParameters.txt"), 
                                      header = TRUE, sep = ",") %>%
    popparam_cleaning() %>%
    mutate(sim_rep = j) %>%
    rename(mean_est2 = value)
  
  point_est_initE_list[[j]] <- point_est_initE
  
  point_est_short_initE <- read.table(paste0(getwd(), "/SEIR_Simulx_new4_pe_short_initE", j, "/populationParameters.txt"), 
                                header = TRUE, sep = ",") %>%
    popparam_cleaning() %>%
    mutate(sim_rep = j) %>%
    rename(mean_est2 = value)
  
  point_est_short_initE_list[[j]] <- point_est_short_initE
}

point_est_df <- do.call("rbind.data.frame", point_est_list) %>%
  mutate(true_value = ifelse(parameter == "beta_ld1", -1.45, -0.5), 
         model = "long")
point_est_df_short <- do.call("rbind.data.frame", point_est_short_list) %>%
  mutate(true_value = ifelse(parameter == "beta_ld1", -1.45, -0.5), 
         model = "short")
point_est_df_short_initE <- do.call("rbind.data.frame", point_est_short_initE_list) %>%
  mutate(true_value = ifelse(parameter == "beta_ld1", -1.45, -0.5), 
         model = "short initE")

comp_df_SEIR4_short_pe <- rbind(point_est_df, point_est_df_short, point_est_df_short_initE)

ggplot(comp_df_SEIR4_short_pe, aes(x = sim_rep, y = mean_est2, col = model)) + 
  geom_point() + 
  facet_wrap(~parameter, ncol = 1, scales = "free_y") +
  geom_line(aes(y = true_value, x = sim_rep), linetype = "dashed", col = "black")



#### bootstrap short ####
SEIR_boot_list_all4_short <- list()
point_est_short_list <- list()
popparams <- c()

# load bootstrap population parameters and calculate bootstrap CIs in 2 ways
for(j in 1:15){
  params_sim_rep <- list()
  
  # load population parameters from individual bootstrap runs 
  for(i in 1:100){
    rm(popparams)
    try(popparams <- read.table(paste0(getwd(), "/SEIR_Simulx_new4_2params_short_", j, "boot", i, "/populationParameters.txt"), 
                                header = TRUE, sep = ",") %>%
          bootstrap_cleaning(., boot_rep = i, sim_rep = j))
    
    
    if(exists("popparams")){
      params_sim_rep[[i]] <- popparams
    }
  }
  
  # calculate CIs
  if(!is_empty(params_sim_rep)){
    pop_params_boot <- bootstrap_CI_calc(params_sim_rep)
  }
  
  
  # load point estimates
  point_est_short <- read.table(paste0(getwd(), "/SEIR_Simulx_new4_pe_short_", j, "/populationParameters.txt"), 
                          header = TRUE, sep = ",") %>%
    popparam_cleaning() %>%
    mutate(sim_rep = j) %>%
    rename(mean_est2 = value)
  
  
  SEIR_boot_list_all4_short[[j]] <- pop_params_boot %>%
    ungroup() %>%
    left_join(point_est_short, by = c("parameter", "sim_rep")) %>%
    select(parameter, sim_rep, mean_est1, mean_est2, CI_LL1, CI_LL2, CI_UL1, CI_UL2, sd_est)
  
}

save(SEIR_boot_list_all4_short, file = "SEIR_boot_list_2params_new4_Simulx_short.RData")



# join point estimates and pivot into longer format 
SEIR_boot_df_short_long4 <- do.call("rbind.data.frame", SEIR_boot_list_all4_short) %>%
  pivot_longer(cols = -c(parameter, sim_rep, sd_est), 
               names_to = c(".value", "method"), 
               names_pattern = "(.*)(\\d+)") %>%
  mutate(parameter = factor(parameter, 
                            levels = c("beta_ld1", "beta_BG1"),
                            labels = c("Lockdown", "Barrier gestures")), 
         true_value = ifelse(parameter == "Lockdown", -1.45, -0.5), 
         method = ifelse(method == 1, "Bootstrap SE", "Empirical bootstrap"))


# plot to compare the two CI methods
ggplot(SEIR_boot_df_short_long4, aes(ymin = CI_LL, ymax = CI_UL, x = sim_rep, y = mean_est, col = method)) + 
  geom_pointrange(position = position_dodge(width = 0.3)) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~parameter, ncol = 1, scale = "free_y") +
  geom_line(aes(y = true_value), linetype = "dashed", col = "darkred", linewidth = 0.8) + 
  labs(title= "Bootstrap SEIR Simulx new4 short estimation", 
       x = "simulation dataset", y = "coefficient value") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2")



# plot comparison between short and long and Simulx 2 values
load("~/PhD/COVID_France/SEIR_vs_Rt_sims/boot_sim_2params_all_new2_Simulx_SEIR/SEIR_boot_list_2params_new2_Simulx.RData")

SEIR_boot_df_short_long_comp <- do.call("rbind.data.frame", SEIR_boot_list_all4_short) %>%
  mutate(model = "45d NPI-free period but 29d cut off") %>%
  bind_rows(do.call("rbind.data.frame", SEIR_boot_list_all4) %>% mutate(model = "45d NPI-free period")) %>%
  bind_rows(do.call("rbind.data.frame", SEIR_boot_list_all) %>% mutate(model = "15d NPI-free period")) %>%
  mutate(parameter = factor(parameter, 
                            levels = c("beta_ld1", "beta_BG1"),
                            labels = c("Lockdown", "Barrier gestures")), 
         true_value = ifelse(parameter == "Lockdown", -1.45, -0.5))


ggplot(SEIR_boot_df_short_long_comp %>% filter(sim_rep <= 20), 
       aes(ymin = CI_LL2, ymax = CI_UL2, x = sim_rep, y = mean_est2, col = model)) + 
  geom_pointrange(position = position_dodge(width = 0.5)) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~parameter, ncol = 1, scale = "free_y") +
  geom_line(aes(y = true_value), linetype = "dashed", col = "darkred", linewidth = 0.8) + 
  labs(title= "Bootstrap SEIR Simulx new4 short vs long estimation", 
       x = "simulation dataset", y = "coefficient value") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom")


## difference in R scripts
file2 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/SEIR_Simulx_boot_2params_new2_1.R"
file4 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/SEIR_Simulx_boot_2params_new4_1.R"

diffr(file2, file4)
