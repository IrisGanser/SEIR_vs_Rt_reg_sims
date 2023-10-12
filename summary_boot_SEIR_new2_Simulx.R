library(tidyverse)
library(colorspace)
library(magrittr)
library(kableExtra)


setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/boot_sim_2params_all_new2_Simulx_SEIR")
source("~/PhD/COVID_France/Dropbox_iris_covid/departement/Données_SPF/Data/data_functions.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")

dir <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/sim_2params_regs"
dir2 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/boot_sim_2params_all_new2_Simulx_SEIRAHD"
dir3 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/boot_point_estimates"

# bootstrap replicates
# j = data simulation
# i = bootstrap replicate

SEIR_boot_list_all <- list()
point_est_list <- list()
popparams <- c()

# load bootstrap population parameters and calculate bootstrap CIs in 2 ways
for(j in 1:100){
  params_sim_rep <- list()
  
  # load population parameters from individual bootstrap runs 
  for(i in 1:100){
    rm(popparams)
    try(popparams <- read.table(paste0(getwd(), "/SEIR_Simulx_new2_2params_", j, "boot", i, "/populationParameters.txt"), 
                                header = TRUE, sep = ",") %>%
          bootstrap_cleaning(., boot_rep = i, sim_rep = j))
    
    
    if(exists("popparams")){
      params_sim_rep[[i]] <- popparams
    }
  }
  
  # calculate CIs
  pop_params_boot <- bootstrap_CI_calc(params_sim_rep)
  
  
  # load point estimates
  point_est <- read.table(paste0(dir3, "/SEIR_Simulx_new2_pe_", j, "/populationParameters.txt"), 
                                    header = TRUE, sep = ",") %>%
    popparam_cleaning() %>%
    mutate(sim_rep = j) %>%
    rename(mean_est2 = value)
  
  
  SEIR_boot_list_all[[j]] <- pop_params_boot %>%
    ungroup() %>%
    left_join(point_est, by = c("parameter", "sim_rep")) %>%
    select(parameter, sim_rep, mean_est1, mean_est2, CI_LL1, CI_LL2, CI_UL1, CI_UL2, sd_est)
    
}

save(SEIR_boot_list_all, file = "SEIR_boot_list_2params_new2_Simulx.RData")



# join point estimates and pivot into longer format 
SEIR_boot_df_long <- do.call("rbind.data.frame", SEIR_boot_list_all) %>%
  pivot_longer(cols = -c(parameter, sim_rep, sd_est), 
               names_to = c(".value", "method"), 
               names_pattern = "(.*)(\\d+)") %>%
  mutate(parameter = factor(parameter, 
                            levels = c("beta_ld1", "beta_BG1"),
                            labels = c("Lockdown", "Barrier gestures")), 
         true_value = ifelse(parameter == "Lockdown", -1.45, -0.5), 
         method = ifelse(method == 1, "Bootstrap SE", "Empirical bootstrap"))


# plot to compare the two CI methods
ggplot(SEIR_boot_df_long, aes(ymin = CI_LL, ymax = CI_UL, x = sim_rep, y = mean_est, col = method)) + 
  geom_pointrange(position = position_dodge(width = 0.3)) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~parameter, ncol = 1, scale = "free_y") +
  geom_line(aes(y = true_value), linetype = "dashed", col = "darkred", linewidth = 0.8) + 
  labs(title= "Bootstrap SEIR Simulx new2", 
       x = "simulation dataset", y = "coefficient value") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2")



load("SEIR_boot_list_2params_new2_Simulx.RData")
load(paste0(dir2, "/SEIRAHD_boot_list_2params_new2_Simulx.RData"))

SEIR_boot_df_comp <- bootstrap_summary(SEIR_boot_list_all) 
SEIRAHD_boot_df_comp <- bootstrap_summary(SEIRAHD_boot_list_all) 



# comparison with regressions
load(paste0(dir, "/reg_res_I_2params_all_Simulx_df.RData"))
load(paste0(dir, "/reg_res_H_2params_all_Simulx_df.RData"))

reg_res_2params_Simulx_new2 <- reg_res_I_2params_all_Simulx_df %>%
  mutate(model = "Regression model IncI") %>%
  bind_rows(reg_res_H_2params_all_Simulx_df %>% mutate(model = "Regression model IncH")) %>%
  reg_summary()


comp_df_Simulx <- SEIR_boot_df_comp %>%
  mutate(model = "SEIR model") %>%
  bind_rows(SEIRAHD_boot_df_comp %>% mutate(model = "SEIRAHD model")) %>% 
  rename(mean_est = mean_est2, CI_LL = CI_LL2, CI_UL = CI_UL2) %>%
  bind_rows(reg_res_2params_Simulx_new2 %>% rename(sim_rep = rep, mean_est = value))


# plot
div_palette <- diverging_hcl("Blue-Red", n = 20)
comp_cols <- c(div_palette[c(1, 5)], "black", div_palette[c(16, 20)])

ggplot(comp_df_Simulx, aes(x = sim_rep, y = mean_est, ymin = CI_LL, ymax = CI_UL, col = model)) + 
  geom_pointrange(position = position_dodge(width = 0.8)) +
  geom_line(aes(y = true_value), col = "black", linetype = "dashed", linewidth = 0.8) + 
  facet_wrap(~parameter, scales = "free_y", nrow = 2) +
  labs(title = "95% CIs, SEIR-type models vs. 2-step regression models", 
       x = "Simulation data set", y = "coefficient value", col = "") + 
  scale_color_manual(values = comp_cols[-3]) +
  theme_bw()


# bring into right format for metrics table
metric_df_Simulx <- comp_df_Simulx %>%
  select(parameter, model, perc_CI_covers, mean_bias, mean_rel_bias) %>%
  unique() %>%
  rename(coverage = perc_CI_covers, bias = mean_bias, relbias = mean_rel_bias) %>%
  pivot_longer(cols = c(coverage, bias, relbias), 
               names_to = "metric", 
               values_to = "value") %>%
  pivot_wider(names_from = model, values_from = value) %>%
  mutate(metric = factor(metric, levels = c("coverage", "bias", "relbias"))) %>%
  arrange (metric, parameter) %>%
  mutate(metric = case_when(metric == "bias" ~ "absolute bias", 
                            metric == "coverage" ~ "CI coverage (%)", 
                            metric == "relbias" ~ "relative bias (%)"))

# print metrics table
metric_df_Simulx %>% 
  kable(digits = 2, format = "html") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE)

