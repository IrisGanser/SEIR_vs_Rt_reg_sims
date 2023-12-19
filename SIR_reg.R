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
library(ggside)

setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/Rt_trajectories")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/deSolve_SIR_model_function.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")

load("sim_SIR_res_df.RData")

ggplot(sim_SIR_res_df, aes(x = time, y = IncI, group = dept_id)) + 
  geom_line()

true_Rt_df <- sim_SIR_res_df %>% dplyr::select(dept_id, time, Rt) %>% rename(day = time, Rt_real = Rt)


reg_data_SIR <- sim_SIR_res_df %>%
  select(dept_id, time, IncI, lockdown1, BG1) %>%
  rename(day = time)



set.seed(123)
selected_depts <- sample(1:94, size = 16, replace = FALSE)

### short NPI-free period ####
cl <- makeCluster(8)
registerDoParallel(cl)


Rt_comp_res_SIR <- Rt_calc_fun(Inc_df = reg_data_SIR, id_col_name = "dept_id", time_col_name = "day", 
                           Inc_col_name = "IncI", model_name = "SIR model", 
                           Rt_ref_df = true_Rt_df, 
                           meansi = 5, stdsi = 5, Rt_prior = 1, Rt_sd_prior = 2)

reg_res_SIR <- Rt_reg_only_fun(data_for_est = Rt_comp_res_SIR$Rt_comp %>%
                                 mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
                                        BG1 = ifelse(day > 70, 1, 0)), 
                       model_name = "SIR model", fits = TRUE)


stopCluster(cl)


sim_SIR_res_df_fitted <- sim_SIR_res_df %>%
  rename(Rt_real = Rt) %>%
  select(-c(lockdown1, BG1)) %>%
  left_join(reg_res_SIR$fits %>% select(-Rt_real), by = c("dept_id", "time" = "day")) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94)))
  

ggplot(sim_SIR_res_df_fitted %>% filter(dept_id %in% selected_depts), 
       aes(x = time)) +
  geom_line(aes(y = Rt, linetype = "Epi Estim")) +
  geom_line(aes(y = Rt_real, linetype = "Real")) + 
  geom_line(aes(y = Rt_fitted, linetype = "Regression fit")) + 
  geom_xsideline(aes(y = IncI)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_wrap(~dept_id2) + 
  labs(#title = "Regression test on SIR-simulated data", 
       linetype = expression(R[t]~ "type"), x = "Day", y = expression(R[t])) +
  scale_x_continuous(expand = c(0.05, 0.05)) + 
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) + 
  theme_bw() + 
  theme(legend.text = element_text(family = "serif", size = 12, hjust = 0), 
        legend.title = element_text(family = "serif", size = 13),
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))


ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/SIR_short_example.jpeg", 
       dpi = 400, width = 14, height = 8.5)


reg_res_SIR$coefs %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.5), 
         bias = value - true_value, 
         rel_bias = abs((value - true_value)/true_value)*100) %>%
  select(-c(lag, model))



### BG1=0.8 short NPI-free period ####
load("sim_SIR_res_df_0.8.RData")

true_Rt_df_0.8 <- sim_SIR_res_df %>% dplyr::select(dept_id, time, Rt) %>% rename(day = time, Rt_real = Rt)

reg_data_SIR_0.8 <- sim_SIR_res_df %>%
  select(dept_id, time, IncI, lockdown1, BG1) %>%
  rename(day = time)


cl <- makeCluster(6)
registerDoParallel(cl)


Rt_comp_res_SIR_0.8 <- Rt_calc_fun(Inc_df = reg_data_SIR_0.8, id_col_name = "dept_id", time_col_name = "day", 
                               Inc_col_name = "IncI", model_name = "SIR model", 
                               Rt_ref_df = true_Rt_df_0.8, 
                               meansi = 5, stdsi = 5, Rt_prior = 1, Rt_sd_prior = 2)

reg_res_SIR_0.8 <- Rt_reg_only_fun(data_for_est = Rt_comp_res_SIR_0.8$Rt_comp %>%
                                 mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
                                        BG1 = ifelse(day > 70, 1, 0)), 
                               model_name = "SIR model", fits = TRUE)


stopCluster(cl)


sim_SIR_res_df_fitted_0.8 <- sim_SIR_res_df %>%
  rename(Rt_real = Rt) %>%
  select(-c(lockdown1, BG1)) %>%
  left_join(reg_res_SIR_0.8$fits %>% select(-Rt_real), by = c("dept_id", "time" = "day")) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94)))


ggplot(sim_SIR_res_df_fitted_0.8 %>% filter(dept_id %in% selected_depts), 
       aes(x = time)) +
  geom_line(aes(y = Rt, linetype = "Epi Estim")) +
  geom_line(aes(y = Rt_real, linetype = "Real")) + 
  geom_line(aes(y = Rt_fitted, linetype = "Regression fit")) + 
  geom_xsideline(aes(y = IncI)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_wrap(~dept_id2) + 
  labs(#title = "Regression test on SIR-simulated data", 
    linetype = expression(R[t]~ "type"), x = "Day", y = expression(R[t])) +
  scale_x_continuous(expand = c(0.05, 0.05)) + 
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) + 
  theme_bw() + 
  theme(legend.text = element_text(family = "serif", size = 12, hjust = 0), 
        legend.title = element_text(family = "serif", size = 13),
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))


ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/SIR_0.8_short_example.jpeg", 
       dpi = 400, width = 14, height = 8.5)


reg_res_SIR_0.8$coefs %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.8), 
         bias = value - true_value, 
         rel_bias = abs((value - true_value)/true_value)*100) %>%
  select(-c(lag, model))


### BG1=0.8 shorter LD short NPI-free period ####
load("sim_SIR_res_df_shortLD_0.8.RData")

true_Rt_df_shortLD_0.8 <- sim_SIR_res_df_shortLD_0.8 %>% dplyr::select(dept_id, time, Rt) %>% rename(day = time, Rt_real = Rt)

reg_data_SIR_shortLD_0.8 <- sim_SIR_res_df_shortLD_0.8 %>%
  select(dept_id, time, IncI, lockdown1, BG1) %>%
  rename(day = time)

new_tstart <- 2:121
new_tend <- 2:121

cl <- makeCluster(6)
registerDoParallel(cl)


Rt_comp_res_SIR_shortLD_0.8 <- Rt_calc_fun(Inc_df = reg_data_SIR_shortLD_0.8, 
                                          id_col_name = "dept_id", time_col_name = "day", 
                                   Inc_col_name = "IncI", model_name = "SIR model", 
                                   Rt_ref_df = true_Rt_df_shortLD_0.8, tstart = new_tstart, tend = new_tend,
                                   meansi = 5, stdsi = 5, Rt_prior = 1, Rt_sd_prior = 2)

reg_res_SIR_shortLD_0.8 <- Rt_reg_only_fun(data_for_est = Rt_comp_res_SIR_shortLD_0.8$Rt_comp %>%
                                     mutate(lockdown1 = ifelse(between(day, 16, 46), 1, 0), 
                                            BG1 = ifelse(day > 46, 1, 0)), 
                                   model_name = "SIR model", fits = TRUE)


stopCluster(cl)


sim_SIR_res_df_fitted_shortLD_0.8 <- sim_SIR_res_df_shortLD_0.8 %>%
  rename(Rt_real = Rt) %>%
  select(-c(lockdown1, BG1)) %>%
  left_join(reg_res_SIR_shortLD_0.8$fits %>% select(-Rt_real), by = c("dept_id", "time" = "day")) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94)))


ggplot(sim_SIR_res_df_fitted_shortLD_0.8 %>% filter(dept_id %in% selected_depts), 
       aes(x = time)) +
  geom_line(aes(y = Rt, linetype = "Epi Estim")) +
  geom_line(aes(y = Rt_real, linetype = "Real")) + 
  geom_line(aes(y = Rt_fitted, linetype = "Regression fit")) + 
  geom_xsideline(aes(y = IncI)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_wrap(~dept_id2) + 
  labs(#title = "Regression test on SIR-simulated data", 
    linetype = expression(R[t]~ "type"), x = "Day", y = expression(R[t])) +
  scale_x_continuous(expand = c(0.05, 0.05)) + 
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) + 
  theme_bw() + 
  theme(legend.text = element_text(family = "serif", size = 12, hjust = 0), 
        legend.title = element_text(family = "serif", size = 13),
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))


ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/SIR_shortLD_0.8_short_example.jpeg", 
       dpi = 400, width = 14, height = 8.5)


reg_res_SIR_shortLD_0.8$coefs %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.8), 
         bias = value - true_value, 
         rel_bias = abs((value - true_value)/true_value)*100) %>%
  select(-c(lag, model))


### longer NPI-free period ####
load("sim_SIR_long_res_df.RData")

true_Rt_df_long <- sim_SIR_long_res_df %>% 
  dplyr::select(dept_id, time, Rt) %>% 
  rename(day = time, Rt_real = Rt)


reg_data_SIR_long <- sim_SIR_long_res_df %>%
  select(dept_id, time, IncI, lockdown1, BG1) %>%
  rename(day = time)

cl <- makeCluster(8)
registerDoParallel(cl)


Rt_comp_res_SIR_long <- Rt_calc_fun(Inc_df = reg_data_SIR_long, id_col_name = "dept_id", time_col_name = "day", 
                               Inc_col_name = "IncI", model_name = "SIR model long", 
                               Rt_ref_df = true_Rt_df_long, 
                               meansi = 5, stdsi = 5, Rt_prior = 1, Rt_sd_prior = 2)

reg_res_SIR_long <- Rt_reg_only_fun(data_for_est = Rt_comp_res_SIR_long$Rt_comp %>%
                                 mutate(lockdown1 = ifelse(between(day, 45, 99), 1, 0), 
                                        BG1 = ifelse(day > 99, 1, 0)), 
                               model_name = "SIR model long", lag_NPIs = FALSE, fits = TRUE)

stopCluster(cl)


sim_SIR_res_df_long_fitted <- sim_SIR_long_res_df %>%
  rename(Rt_real = Rt) %>%
  select(-c(lockdown1, BG1)) %>%
  left_join(reg_res_SIR_long$fits %>% select(-Rt_real), by = c("dept_id", "time" = "day")) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94)))


ggplot(sim_SIR_res_df_long_fitted %>% filter(dept_id %in% selected_depts), 
       aes(x = time)) +
  geom_line(aes(y = Rt, linetype = "Epi Estim")) +
  geom_line(aes(y = Rt_real, linetype = "Real")) + 
  geom_line(aes(y = Rt_fitted, linetype = "Regression")) + 
  geom_xsideline(aes(y = IncI)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_wrap(~dept_id2) +
  labs(#title = "Regression test on SIR-simulated data (longer NPI-free period)", 
  linetype = expression(R[t]~ "type"), x = "Day", y = expression(R[t])) +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) + 
  theme_bw() + 
  theme(legend.text = element_text(family = "serif", size = 12, hjust = 0), 
        legend.title = element_text(family = "serif", size = 13),
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))


ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/SIR_long_example.jpeg", 
       dpi = 400, width = 14, height = 8.5)

reg_res_SIR_long$coefs %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.5), 
         bias = value - true_value, 
         rel_bias = abs((value - true_value)/true_value)*100) %>%
  select(-c(lag, model))




### BG1=0.8 longer NPI-free period ####
load("sim_SIR_long_res_df_0.8.RData")

true_Rt_df_long_0.8 <- sim_SIR_long_res_df %>% 
  dplyr::select(dept_id, time, Rt) %>% 
  rename(day = time, Rt_real = Rt)


reg_data_SIR_long_0.8 <- sim_SIR_long_res_df %>%
  select(dept_id, time, IncI, lockdown1, BG1) %>%
  rename(day = time)

cl <- makeCluster(8)
registerDoParallel(cl)


Rt_comp_res_SIR_long_0.8 <- Rt_calc_fun(Inc_df = reg_data_SIR_long_0.8, id_col_name = "dept_id", time_col_name = "day", 
                                    Inc_col_name = "IncI", model_name = "SIR model long", 
                                    Rt_ref_df = true_Rt_df_long_0.8, 
                                    meansi = 5, stdsi = 5, Rt_prior = 1, Rt_sd_prior = 2)

reg_res_SIR_long_0.8 <- Rt_reg_only_fun(data_for_est = Rt_comp_res_SIR_long_0.8$Rt_comp %>%
                                      mutate(lockdown1 = ifelse(between(day, 45, 99), 1, 0), 
                                             BG1 = ifelse(day > 99, 1, 0)), 
                                    model_name = "SIR model long", lag_NPIs = FALSE, fits = TRUE)

stopCluster(cl)


sim_SIR_res_df_long_fitted_0.8 <- sim_SIR_long_res_df %>%
  rename(Rt_real = Rt) %>%
  select(-c(lockdown1, BG1)) %>%
  left_join(reg_res_SIR_long_0.8$fits %>% select(-Rt_real), by = c("dept_id", "time" = "day")) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94)))


ggplot(sim_SIR_res_df_long_fitted_0.8 %>% filter(dept_id %in% selected_depts), 
       aes(x = time)) +
  geom_line(aes(y = Rt, linetype = "Epi Estim")) +
  geom_line(aes(y = Rt_real, linetype = "Real")) + 
  geom_line(aes(y = Rt_fitted, linetype = "Regression")) + 
  geom_xsideline(aes(y = IncI)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_wrap(~dept_id2) +
  labs(#title = "Regression test on SIR-simulated data (longer NPI-free period)", 
    linetype = expression(R[t]~ "type"), x = "Day", y = expression(R[t])) +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) + 
  theme_bw() + 
  theme(legend.text = element_text(family = "serif", size = 12, hjust = 0), 
        legend.title = element_text(family = "serif", size = 13),
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))


ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/SIR_0.8_long_example.jpeg", 
       dpi = 400, width = 14, height = 8.5)

reg_res_SIR_long_0.8$coefs %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.8), 
         bias = value - true_value, 
         rel_bias = abs((value - true_value)/true_value)*100) %>%
  select(-c(lag, model))


# shorter LD1
new_tstart <- 2:121
new_tend <- 2:121

load("sim_SIR_long_res_df_shortLD_0.8.RData")

true_Rt_df_long_shortLD_0.8 <- sim_SIR_long_res_df_shortLD_0.8 %>% 
  filter(time <= 121) %>%
  dplyr::select(dept_id, time, Rt) %>% 
  rename(day = time, Rt_real = Rt)


reg_data_SIR_long_shortLD_0.8 <- sim_SIR_long_res_df_shortLD_0.8 %>%
  filter(time <= 121) %>%
  select(dept_id, time, IncI, lockdown1, BG1) %>%
  rename(day = time)

cl <- makeCluster(8)
registerDoParallel(cl)


Rt_comp_res_SIR_long_shortLD_0.8 <- Rt_calc_fun(Inc_df = reg_data_SIR_long_shortLD_0.8, 
                                                id_col_name = "dept_id", time_col_name = "day", 
                                                Inc_col_name = "IncI", model_name = "SIR model long", 
                                                Rt_ref_df = true_Rt_df_long_shortLD_0.8, 
                                                tstart = new_tstart, tend = new_tend,
                                                meansi = 5, stdsi = 5, Rt_prior = 1, Rt_sd_prior = 2)

reg_res_SIR_long_shortLD_0.8 <- Rt_reg_only_fun(data_for_est = Rt_comp_res_SIR_long_shortLD_0.8$Rt_comp %>%
                                          mutate(lockdown1 = ifelse(between(day, 45, 65), 1, 0), 
                                                 BG1 = ifelse(day > 85, 1, 0)), 
                                        model_name = "SIR model long", lag_NPIs = FALSE, fits = TRUE)

stopCluster(cl)


sim_SIR_res_df_long_fitted_shortLD_0.8 <- sim_SIR_long_res_df_shortLD_0.8 %>%
  filter(time <= 121) %>%
  rename(Rt_real = Rt) %>%
  select(-c(lockdown1, BG1)) %>%
  left_join(reg_res_SIR_long_shortLD_0.8$fits %>% select(-Rt_real), by = c("dept_id", "time" = "day")) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94)))


ggplot(sim_SIR_res_df_long_fitted_shortLD_0.8 %>% filter(dept_id %in% selected_depts), 
       aes(x = time)) +
  geom_line(aes(y = Rt, linetype = "Epi Estim")) +
  geom_line(aes(y = Rt_real, linetype = "Real")) + 
  geom_line(aes(y = Rt_fitted, linetype = "Regression")) + 
  geom_xsideline(aes(y = IncI)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_wrap(~dept_id2) +
  labs(#title = "Regression test on SIR-simulated data (longer NPI-free period)", 
    linetype = expression(R[t]~ "type"), x = "Day", y = expression(R[t])) +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) + 
  theme_bw() + 
  theme(legend.text = element_text(family = "serif", size = 12, hjust = 0), 
        legend.title = element_text(family = "serif", size = 13),
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))


ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/SIR_shortLD_0.8_long_example.jpeg", 
       dpi = 400, width = 14, height = 8.5)

reg_res_SIR_long_shortLD_0.8$coefs %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.8), 
         bias = value - true_value, 
         rel_bias = abs((value - true_value)/true_value)*100) %>%
  select(-c(lag, model))




### longer NPI-free period m10 ####
reg_res_SIR_long_m10 <- Rt_reg_only_fun(data_for_est = Rt_comp_res_SIR_long$Rt_comp %>%
                                          mutate(lockdown1 = ifelse(between(day, 45, 99), 1, 0), 
                                                 BG1 = ifelse(day > 99, 1, 0)) %>%
                                          filter(day > 13), 
                                        model_name = "SIR model long m10", 
                                        lag_NPIs = FALSE, fits = TRUE)



sim_SIR_res_df_long_m10_fitted <- sim_SIR_long_res_df %>%
  rename(Rt_real = Rt) %>%
  select(-c(lockdown1, BG1)) %>%
  left_join(reg_res_SIR_long_m10$fits %>% select(-Rt_real), by = c("dept_id", "time" = "day")) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94)))

ggplot(sim_SIR_res_df_long_m10_fitted %>% filter(dept_id %in% selected_depts), 
       aes(x = time)) +
  geom_line(aes(y = Rt, linetype = "Epi Estim")) +
  geom_line(aes(y = Rt_real, linetype = "Real")) + 
  geom_line(aes(y = Rt_fitted, linetype = "Regression")) + 
  geom_xsideline(aes(y = IncI)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_wrap(~dept_id2) +
  labs(linetype = "Rt type", title = "Regression test on SIR-simulated data (longer NPI-free period), 10 days cut off") +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) + 
  theme_bw() +
  theme(legend.text = element_text(family = "serif", size = 12, hjust = 0), 
        legend.title = element_blank(),
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))

reg_res_SIR_long_m10$coefs %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.5), 
         bias = value - true_value, 
         rel_bias = abs((value - true_value)/true_value)*100) %>%
  select(-c(lag, model))




### BG1=0.8 longer NPI-free period m10 ####
reg_res_SIR_long_m10_0.8 <- Rt_reg_only_fun(data_for_est = Rt_comp_res_SIR_long_0.8$Rt_comp %>%
                                          mutate(lockdown1 = ifelse(between(day, 45, 99), 1, 0), 
                                                 BG1 = ifelse(day > 99, 1, 0)) %>%
                                          filter(day > 13), 
                                        model_name = "SIR model long m10", 
                                        lag_NPIs = FALSE, fits = TRUE)



sim_SIR_res_df_long_m10_fitted_0.8 <- sim_SIR_long_res_df %>%
  rename(Rt_real = Rt) %>%
  select(-c(lockdown1, BG1)) %>%
  left_join(reg_res_SIR_long_m10_0.8$fits %>% select(-Rt_real), by = c("dept_id", "time" = "day")) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94)))

ggplot(sim_SIR_res_df_long_m10_fitted_0.8 %>% filter(dept_id %in% selected_depts), 
       aes(x = time)) +
  geom_line(aes(y = Rt, linetype = "Epi Estim")) +
  geom_line(aes(y = Rt_real, linetype = "Real")) + 
  geom_line(aes(y = Rt_fitted, linetype = "Regression")) + 
  geom_xsideline(aes(y = IncI)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_wrap(~dept_id2) +
  labs(linetype = "Rt type", title = "Regression test on SIR-simulated data (longer NPI-free period), 10 days cut off") +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) + 
  theme_bw() +
  theme(legend.text = element_text(family = "serif", size = 12, hjust = 0), 
        legend.title = element_blank(),
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))

reg_res_SIR_long_m10_0.8$coefs %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.8), 
         bias = value - true_value, 
         rel_bias = abs((value - true_value)/true_value)*100) %>%
  select(-c(lag, model))



# shorter LD1
reg_res_SIR_long_m10_shortLD_0.8 <- Rt_reg_only_fun(data_for_est = Rt_comp_res_SIR_long_shortLD_0.8$Rt_comp %>%
                                              mutate(lockdown1 = ifelse(between(day, 45, 65), 1, 0), 
                                                     BG1 = ifelse(day > 85, 1, 0)) %>%
                                              filter(day > 11), 
                                            model_name = "SIR model long m10", 
                                            lag_NPIs = FALSE, fits = TRUE)



sim_SIR_res_df_long_m10_fitted_shortLD_0.8 <- sim_SIR_long_res_df_shortLD_0.8 %>%
  filter(time <=121) %>%
  rename(Rt_real = Rt) %>%
  select(-c(lockdown1, BG1)) %>%
  left_join(reg_res_SIR_long_m10_shortLD_0.8$fits %>% select(-Rt_real), by = c("dept_id", "time" = "day")) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94)))

ggplot(sim_SIR_res_df_long_m10_fitted_shortLD_0.8 %>% filter(dept_id %in% selected_depts), 
       aes(x = time)) +
  geom_line(aes(y = Rt, linetype = "Epi Estim")) +
  geom_line(aes(y = Rt_real, linetype = "Real")) + 
  geom_line(aes(y = Rt_fitted, linetype = "Regression")) + 
  geom_xsideline(aes(y = IncI)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_wrap(~dept_id2) +
  labs(linetype = "Rt type", title = "Regression test on SIR-simulated data (longer NPI-free period), 10 days cut off") +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) + 
  theme_bw() +
  theme(legend.text = element_text(family = "serif", size = 12, hjust = 0), 
        legend.title = element_blank(),
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))

ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/SIR_m10_shortLD_0.8_long_example.jpeg", 
       dpi = 400, width = 14, height = 8.5)

reg_res_SIR_long_m10_shortLD_0.8$coefs %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.8), 
         bias = value - true_value, 
         rel_bias = abs((value - true_value)/true_value)*100) %>%
  select(-c(lag, model))




# Metrics table -----------------------------------------------------------

metrics_table <- reg_res_SIR_long_m10_shortLD_0.8$coefs %>%
  mutate(model = "late NPI 1, 10 days cut off") %>%
  bind_rows(reg_res_SIR_long_shortLD_0.8$coefs %>% mutate(model = "late NPI 1")) %>%
  bind_rows(reg_res_SIR_shortLD_0.8$coefs %>% mutate(model = "early NPI 1")) %>%
  mutate(reg_coef = paste0(round(value, 3), " [", 
                           round(CI_LL, 3), ", ", round(CI_UL, 3), "]")) %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.8), 
         bias = as.character(round(value - true_value, 3)), 
         relbias = as.character(round(abs((value - true_value)/true_value)*100, 3))) %>%
  select(-c(lag, value, CI_LL, CI_UL, true_value)) %>%
  pivot_longer(cols = c(reg_coef, bias, relbias), 
             names_to = "metric", 
             values_to = "value") %>%
  pivot_wider(names_from = model, values_from = value, 
              names_sort = TRUE) %>%
  mutate(metric = factor(metric, levels = c("reg_coef", "bias", "relbias"))) %>%
  arrange (parameter, metric) %>%
  mutate(metric = case_when(metric == "bias" ~ "absolute bias", 
                            metric == "reg_coef" ~ "regression coefficient", 
                            metric == "relbias" ~ "relative bias (%)"))

# print metrics table
metrics_table %>% 
  ungroup() %>%
  dplyr::select(-parameter) %>%
  kable(digits = 2, format = "latex" #, table.attr = "style='width:40%;'"
  ) %>%
  kable_styling(bootstrap_options = "striped") %>%
  pack_rows("NPI 1", 1, 3) %>%
  pack_rows("NPI 2", 4, 6) %>%
  column_spec(1, width_min = "3.8cm")

