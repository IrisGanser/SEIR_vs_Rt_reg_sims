library(tidyverse)
library(deSolve)
library(lme4)
library(EpiEstim)
library(colorspace)
library(parallel)
library(foreach)
library(doParallel)
library(magrittr)
library(RColorBrewer)


setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/Rt_trajectories")
source("~/PhD/COVID_France/Dropbox_iris_covid/departement/Données_SPF/Data/data_functions.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")
dir1 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"


# load necessary data
load("sim_res_ld1_start.RData")
load("sim_res_ld1_start_high.RData")
load("sim_res_ld1_strength.RData")

popsize_df <- read.csv(paste(dir1, "popsize_df.csv", sep = "/")) %>%
  mutate(dept_id = ifelse(dept_id < 20, dept_id, dept_id - 1))

ind_params_lowb1 <- read.table("ind_params_lowb1.txt", header = TRUE, sep = " ")
ind_params <- read.table("ind_params.txt", header = TRUE, sep = " ")


ld1_start <- c(20, 30, 40, 50, 60)
ld1_strength <- seq(0.5, 2, 0.1)

# EpiEstim
cl <- makeCluster(10)
registerDoParallel(cl)

res_EE_ld1_start_list <- list()

for(j in 1:length(sim_res_ld1_start)){
  data_x <- sim_res_ld1_start[[j]] %>%
    left_join(popsize_df, by = c("id" = "dept_id")) %>%
    left_join(ind_params_lowb1, by = "id") %>%
    mutate(Rt_SEIRAHD = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
           IncI_unscaled = round(IncI*popsize/10^4), 
           lockdown1 = ifelse(between(time, ld1_start[j], 70), 1, 0), 
           BG1 = ifelse(time > 70, 1, 0)) %>% 
    rename(dept_id = id, day = time)
  
  res_EE_ld1_start <- EpiEstim_only_fun(data_for_est = data_x, 
                                        Inc_name = "IncI_unscaled", 
                                        meansi = 10.1, stdsi = 8.75) %>%
    mutate(ld1_start = ld1_start[j]) %>%
    left_join(data_x, by = c("dept_id", "day"))
  
  res_EE_ld1_start_list[[j]] <- res_EE_ld1_start
}



res_EE_ld1_start_high_list <- list()

for(j in 1:length(sim_res_ld1_start)){
  data_x <- sim_res_ld1_start_high[[j]] %>%
    left_join(popsize_df, by = c("id" = "dept_id")) %>%
    left_join(ind_params, by = "id") %>%
    mutate(Rt_SEIRAHD = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
           IncI_unscaled = round(IncI*popsize/10^4), 
           lockdown1 = ifelse(between(time, ld1_start[j], 70), 1, 0), 
           BG1 = ifelse(time > 70, 1, 0)) %>% 
    rename(dept_id = id, day = time)
  
  res_EE_ld1_start_high <- EpiEstim_only_fun(data_for_est = data_x, 
                                             Inc_name = "IncI_unscaled", 
                                             meansi = 10.1, stdsi = 8.75) %>%
    mutate(ld1_start = ld1_start[j]) %>%
    left_join(data_x, by = c("dept_id", "day"))
  
  res_EE_ld1_start_high_list[[j]] <- res_EE_ld1_start_high
}


res_EE_ld1_strength_list <- list()

for(j in 1:length(sim_res_ld1_strength)){
  ind_prms <- read.table(paste0("ind_params", ld1_strength[j], ".txt"), header = TRUE, sep = " ")
  data_x <- sim_res_ld1_strength[[j]] %>%
    left_join(popsize_df, by = c("id" = "dept_id")) %>%
    left_join(ind_prms, by = "id") %>%
    mutate(Rt_SEIRAHD = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
           IncI_unscaled = round(IncI*popsize/10^4), 
           lockdown1 = ifelse(between(time, 20, 70), 1, 0), 
           BG1 = ifelse(time > 70, 1, 0)) %>% 
    rename(dept_id = id, day = time)
  
  res_EE_ld1_strength <- EpiEstim_only_fun(data_for_est = data_x, 
                                           Inc_name = "IncI_unscaled", 
                                           meansi = 10.1, stdsi = 8.75) %>%
    mutate(ld1_strength = ld1_strength[j]) %>%
    left_join(data_x, by = c("dept_id", "day"))
  
  res_EE_ld1_strength_list[[j]] <- res_EE_ld1_strength
}

stopCluster(cl)


# compare and plot results
res_EE_ld1_start_df <- do.call("rbind.data.frame", res_EE_ld1_start_list) 

res_EE_ld1_start_high_df <- do.call("rbind.data.frame", res_EE_ld1_start_high_list) 

res_EE_ld1_strength_df <- do.call("rbind.data.frame", res_EE_ld1_strength_list) 


ggplot(res_EE_ld1_start_df %>% filter(dept_id %in% c(1:6)),
       aes(x = day, y = Rt, col = as.factor(ld1_start))) + 
  geom_line() +
  geom_ribbon(aes(ymax = CI_UL, ymin = CI_LL, fill = as.factor(ld1_start)), alpha = 0.5) +
  geom_line(aes(y = Rt_SEIRAHD), linetype = "dashed") + 
  facet_wrap(~dept_id) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  labs(title = "Comparison Rt estimated by EpiEstim and real underlying Rt (lower basic transmission)", 
       col = "Lockdown onset day", fill = "Lockdown onset day") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + 
  theme(legend.position = "bottom")

ggplot(res_EE_ld1_start_high_df %>% filter(dept_id %in% c(1:6)),
       aes(x = day, y = Rt, col = as.factor(ld1_start))) + 
  geom_line() +
  geom_ribbon(aes(ymax = CI_UL, ymin = CI_LL, fill = as.factor(ld1_start)), alpha = 0.5) +
  geom_line(aes(y = Rt_SEIRAHD), linetype = "dashed") + 
  facet_wrap(~dept_id) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  labs(title = "Comparison Rt estimated by EpiEstim and real underlying Rt (higher basic transmission)", 
       col = "Lockdown onset day", fill = "Lockdown onset day") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + 
  theme(legend.position = "bottom")


q4 <- qualitative_hcl(16, palette = "Dark 3")
ggplot(res_EE_ld1_strength_df %>% filter(dept_id == 1),
       aes(x = day, y = Rt, col = as.factor(ld1_strength))) + 
  geom_line(aes(y = Rt_SEIRAHD), linetype = "dashed") + 
  geom_line() +
  geom_ribbon(aes(ymax = CI_UL, ymin = CI_LL, fill = as.factor(ld1_strength)), alpha = 0.5) +
  facet_wrap(~ld1_strength) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  labs(title = "Comparison Rt estimated by EpiEstim and real underlying Rt", 
       col = "Lockdown strength", fill = "Lockdown strength") +
  theme_bw() +
  scale_color_manual(values = q4) +
  scale_fill_manual(values = q4) + 
  theme(legend.position = "bottom")

#regressions
cl <- makeCluster(10)
registerDoParallel(cl)

res_reg_ld1_start_list <- list()

for(j in 1:length(sim_res_ld1_start)){
  data_x <- sim_res_ld1_start[[j]] %>%
    left_join(popsize_df, by = c("id" = "dept_id")) %>%
    left_join(ind_params_lowb1, by = "id") %>%
    mutate(Rt_SEIRAHD = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
           IncI_unscaled = round(IncI*popsize/10^4), 
           lockdown1 = ifelse(between(time, ld1_start[j], 70), 1, 0), 
           BG1 = ifelse(time > 70, 1, 0)) %>% 
    rename(dept_id = id, day = time)
  
  res_reg_ld1_start <- EpiEstim_reg_fun(data_for_est = data_x, 
                                        Inc_name = "IncI_unscaled", 
                                        rep_num = 1, 
                                        meansi = 10.1, stdsi = 8.75) %>%
    mutate(ld1_start = ld1_start[j]) 
  
  res_reg_ld1_start_list[[j]] <- res_reg_ld1_start
}


res_reg_ld1_start_high_list <- list()

for(j in 1:length(sim_res_ld1_start)){
  data_x <- sim_res_ld1_start_high[[j]] %>%
    left_join(popsize_df, by = c("id" = "dept_id")) %>%
    left_join(ind_params, by = "id") %>%
    mutate(Rt_SEIRAHD = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
           IncI_unscaled = round(IncI*popsize/10^4), 
           lockdown1 = ifelse(between(time, ld1_start[j], 70), 1, 0), 
           BG1 = ifelse(time > 70, 1, 0)) %>% 
    rename(dept_id = id, day = time)
  
  res_reg_ld1_start_high <- EpiEstim_reg_fun(data_for_est = data_x, 
                                             Inc_name = "IncI_unscaled", 
                                             rep_num = 1, 
                                             meansi = 10.1, stdsi = 8.75) %>%
    mutate(ld1_start = ld1_start[j]) 
  
  res_reg_ld1_start_high_list[[j]] <- res_reg_ld1_start_high
}


res_reg_ld1_strength_list <- list()

for(j in 1:length(sim_res_ld1_strength)){
  ind_prms <- read.table(paste0("ind_params", ld1_strength[j], ".txt"), header = TRUE, sep = " ")
  data_x <- sim_res_ld1_strength[[j]] %>%
    left_join(popsize_df, by = c("id" = "dept_id")) %>%
    left_join(ind_prms, by = "id") %>%
    mutate(Rt_SEIRAHD = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
           IncI_unscaled = round(IncI*popsize/10^4), 
           lockdown1 = ifelse(between(time, 20, 70), 1, 0), 
           BG1 = ifelse(time > 70, 1, 0)) %>% 
    rename(dept_id = id, day = time)
  
  res_reg_ld1_strength <- EpiEstim_reg_fun(data_for_est = data_x, 
                                           Inc_name = "IncI_unscaled", 
                                           rep_num = 1, 
                                           meansi = 10.1, stdsi = 8.75) %>%
    mutate(ld1_strength = ld1_strength[j]) 
  
  res_reg_ld1_strength_list[[j]] <- res_reg_ld1_strength
}

stopCluster(cl)


res_reg_ld1_start_df <- do.call("rbind.data.frame", res_reg_ld1_start_list) %>%
  mutate(true_value = ifelse(parameter == "Lockdown 1", -1.45, -0.5), 
         bias = abs(true_value - value), 
         rel_bias = abs(true_value - value)/abs(true_value)*100)
res_reg_ld1_start_high_df <- do.call("rbind.data.frame", res_reg_ld1_start_high_list) %>%
  mutate(true_value = ifelse(parameter == "Lockdown 1", -1.45, -0.5), 
         bias = abs(true_value - value), 
         rel_bias = abs(true_value - value)/abs(true_value)*100)
res_reg_ld1_strength_df <- do.call("rbind.data.frame", res_reg_ld1_strength_list) %>%
  mutate(true_value = ifelse(parameter == "Lockdown 1", 0-ld1_strength, -0.5), 
         bias = abs(true_value - value), 
         rel_bias = abs(true_value - value)/abs(true_value)*100)


# plots
ggplot(res_reg_ld1_start_df, aes(ymin = CI_LL, ymax = CI_UL, x = ld1_start, 
                                 y = value, col = as.factor(ld1_start))) + 
  geom_pointrange(position = position_dodge(width = 1)) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~parameter, ncol = 1, scale = "free_y") +
  geom_line(aes(y = true_value), linetype = "dashed", col = "black", linewidth = 0.8) + 
  labs(title = "Comparison 2 step regression: NPI start day (low transmission)", col = "Lockdown1 start day",
       x = "Lockdown1 start day", y = "coefficient value") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2")

ggplot(res_reg_ld1_start_high_df, aes(ymin = CI_LL, ymax = CI_UL, x = ld1_start, 
                                      y = value, col = as.factor(ld1_start))) + 
  geom_pointrange(position = position_dodge(width = 1)) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~parameter, ncol = 1, scale = "free_y") +
  geom_line(aes(y = true_value), linetype = "dashed", col = "black", linewidth = 0.8) + 
  labs(title = "Comparison 2 step regression: NPI start day (high transmission)", col = "Lockdown1 start day",
       x = "Lockdown1 start day", y = "coefficient value") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2")


ggplot(res_reg_ld1_strength_df, aes(ymin = CI_LL, ymax = CI_UL, x = ld1_strength, 
                                      y = value, col = as.factor(ld1_strength))) + 
  geom_pointrange() + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~parameter, ncol = 1, scale = "free_y") +
  geom_segment(aes(y = true_value, yend = true_value, x = ld1_strength-0.05, xend = ld1_strength+0.05), 
               linetype = "dashed", col = "black", linewidth = 0.8) + 
  labs(title = "Comparison 2 step regression: NPI strength", col = "Lockdown1 strength",
       x = "Lockdown1 strength", y = "coefficient value") +
  theme_bw() +
  scale_color_manual(values = q4)


# tables
metric_table_reg_ld1_strength <- res_reg_ld1_strength_df %>%
  mutate(parameter = ifelse(parameter == "Lockdown 1", "NPI 1", "NPI 2")) %>%
  select(parameter, ld1_strength, bias, rel_bias) %>%
  unique() %>%
  pivot_longer(cols = c(bias, rel_bias), 
               names_to = "metric", 
               values_to = "value") %>%
  pivot_wider(names_from = ld1_strength, values_from = value) %>%
  mutate(metric = factor(metric, levels = c("bias", "rel_bias"))) %>%
  arrange (parameter, metric) %>%
  mutate(metric = case_when(metric == "bias" ~ "absolute bias", 
                            metric == "rel_bias" ~ "relative bias (%)"))

metric_table_reg_ld1_strength %>% 
  select(-parameter) %>%
  kable(digits = 2, format = "html") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
  pack_rows("NPI 1", 1, 2) %>%
  pack_rows("NPI 2", 3, 4)


# save results
save(res_reg_ld1_start_df, res_reg_ld1_start_high_df, res_reg_ld1_strength_df, 
     res_EE_ld1_start_df , res_EE_ld1_start_high_df, res_EE_ld1_strength_df, 
     file = "Rt_NPI_scen_res_dfs.RData")
