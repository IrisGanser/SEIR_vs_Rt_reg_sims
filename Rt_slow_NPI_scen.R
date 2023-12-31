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



# Simulx 2 ----------------------------------------------------------------

# load necessary data
load("sim_res_slow_NPI.RData")
sim_res_slow_names <- c("1week_slow", "2week_slow", "1week_slow_late", "2week_slow_late")


popsize_df <- read.csv(paste(dir1, "popsize_df.csv", sep = "/")) %>%
  mutate(dept_id = ifelse(dept_id < 20, dept_id, dept_id - 1))

ind_params_lowb1 <- read.table("ind_params_lowb1.txt", header = TRUE, sep = " ")
ind_params <- read.table("ind_params.txt", header = TRUE, sep = " ")


# EpiEstim
cl <- makeCluster(10)
registerDoParallel(cl)

res_EE_slow_NPI_list <- list()
res_reg_slow_NPI_list <- list()
fits_reg_slow_NPI_list <- list()

for(j in 1:4){
  
  if(j %in% 1:2){
    indprms <- ind_params
  }else{
    indprms <- ind_params_lowb1
  }
 
  # assemble data
  data_SEIR_monolix <- read.table(paste0("data_sim_SEIR_Simulx_2params_", sim_res_slow_names[j], ".txt"), 
                                  sep = ",", header = TRUE) %>%
    select(dept_id, day, lockdown1, BG1)
  
  data_x <- sim_res_slow_NPI[[j]] %>%
    left_join(popsize_df, by = c("id" = "dept_id")) %>%
    left_join(indprms, by = "id") %>%
    mutate(Rt_SEIRAHD = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
           IncI_unscaled = round(IncI*popsize/10^4)) %>% 
    rename(dept_id = id, day = time) %>%
    left_join(data_SEIR_monolix, by = c("dept_id", "day"))
  
  # EpiEstim only
  res_EE <- EpiEstim_only_fun(data_for_est = data_x, 
                              Inc_name = "IncI_unscaled", 
                              meansi = 10.1, stdsi = 8.75) %>%
    mutate(model = sim_res_slow_names[j]) %>%
    left_join(data_x, by = c("dept_id", "day"))
  
  res_EE_slow_NPI_list[[j]] <- res_EE
  
  # regression
  reg_res <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = res_EE)
  
  # coefficients
  coefs_Rt_reg <- coefficients(reg_res)$dept_id
  
  confint_Rt_reg <- data.frame(confint(reg_res, method="Wald"))[-c(1:2), ]
  names(confint_Rt_reg) <- c("CI_LL", "CI_UL")
  
  
  coefs_Rt_reg <- coefs_Rt_reg %>%
    dplyr::select(-1) %>%
    unique() %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
    cbind(confint_Rt_reg[-1,]) %>%
    mutate(parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("Lockdown 1", "Barrier gestures")), 
           model = sim_res_slow_names[j])
  
  res_reg_slow_NPI_list[[j]] <- coefs_Rt_reg
  
  
  # regression fits
  fitted_vals <- exp(fitted(reg_res))
  fitted_df <- res_EE %>%
    filter(!is.na(Rt)) %>%
    mutate(Rt_fitted = fitted_vals) %>%
    select(dept_id, day, Rt_SEIRAHD, Rt, CI_LL, CI_UL, Rt_fitted, model, lockdown1, BG1)
  
  fits_reg_slow_NPI_list[[j]] <- fitted_df
}

stopCluster(cl)


# EpiEstim plots
res_EE_slow_NPI_df <- do.call("rbind.data.frame", res_EE_slow_NPI_list) 

ggplot(res_EE_slow_NPI_df %>% filter(dept_id %in% c(1, 4, 9, 13)),
       aes(x = day, y = Rt, col = model)) + 
  geom_line() +
  geom_ribbon(aes(ymax = CI_UL, ymin = CI_LL, fill = model), alpha = 0.5) +
  geom_line(aes(y = Rt_SEIRAHD), linetype = "dashed") + 
  facet_wrap(~dept_id) + 
  scale_x_continuous(breaks = seq(0, 120, 20)) + 
  labs(title = "Comparison Rt estimated by EpiEstim and real underlying Rt", 
       col = "", fill = "") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + 
  theme(legend.position = "bottom")


# regression plots
res_reg_slow_NPI_df <- do.call("rbind.data.frame", res_reg_slow_NPI_list) %>%
  mutate(true_value = ifelse(parameter == "Lockdown 1", -1.45, -0.5), 
         bias = abs(true_value - value), 
         rel_bias = abs(true_value - value)/abs(true_value)*100)

ggplot(res_reg_slow_NPI_df, aes(ymin = CI_LL, ymax = CI_UL, x = model, 
                                 y = value, col = model)) + 
  geom_pointrange(position = position_dodge(width = 1)) +  
  facet_wrap(~parameter, ncol = 1, scale = "free_y") +
  geom_hline(aes(yintercept = true_value), linetype = "dashed", col = "black", linewidth = 0.8) + 
  labs(title = "Comparison 2 step regression: slower NPI onset", col = "",
       x = "", y = "coefficient value") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())


# regression fits
fits_reg_slow_NPI_df <- do.call("rbind.data.frame", fits_reg_slow_NPI_list)

ggplot(fits_reg_slow_NPI_df %>% filter(dept_id %in% c(1, 4, 9, 13)),
       aes(x = day, y = Rt, col = model)) + 
  geom_line(aes(linetype = "EpiEstim Rt"), linewidth = 0.7) +
  geom_ribbon(aes(ymax = CI_UL, ymin = CI_LL, fill = model), alpha = 0.4) +
  geom_line(aes(y = Rt_SEIRAHD, linetype = "True Rt"), linewidth = 0.7) + 
  geom_line(aes(y = Rt_fitted, linetype = "Regression fit Rt"), linewidth = 0.7) + 
  facet_wrap(~dept_id) + 
  scale_x_continuous(breaks = seq(0, 120, 20)) + 
  scale_linetype_manual(values = c("True Rt" = "dashed", 
                                   "EpiEstim Rt" = "solid", 
                                   "Regression fit Rt" = "dotted")) + 
  labs(title = expression(R[t] ~ "estimated by EpiEstim and real underlying" ~ R[t]), 
       col = "", fill = "", linetype = "") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + 
  theme(legend.position = "bottom")



# Simulx 3 ----------------------------------------------------------------

# load necessary data
load("sim_res3_slow_NPI.RData")
sim_res_slow_names <- c("1week_slow", "2week_slow", "1week_slow_late", "2week_slow_late")


popsize_df <- read.csv(paste(dir1, "popsize_df.csv", sep = "/")) %>%
  mutate(dept_id = ifelse(dept_id < 20, dept_id, dept_id - 1))

ind_params3_lowb1 <- read.table("ind_params3_lowb1.txt", header = TRUE, sep = " ")
ind_params3 <- read.table("ind_params3.txt", header = TRUE, sep = " ")


# EpiEstim
cl <- makeCluster(10)
registerDoParallel(cl)

res_EE_slow_NPI_list3 <- list()
res_reg_slow_NPI_list3 <- list()
fits_reg_slow_NPI_list3 <- list()

for(j in 1:4){
  
  # if(j %in% 1:2){
  #   
  # }else{
  #   indprms <- ind_params3_lowb1
  # }
  indprms <- ind_params3
  
  # assemble data
  data_SEIR_monolix <- read.table(paste0("data_sim_SEIR_Simulx3_2params_", sim_res_slow_names[j], ".txt"), 
                                  sep = ",", header = TRUE) %>%
    select(dept_id, day, lockdown1, BG1)
  
  data_x <- sim_res3_slow_NPI[[j]] %>%
    left_join(popsize_df, by = c("id" = "dept_id")) %>%
    left_join(indprms, by = "id") %>%
    mutate(Rt_SEIRAHD = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
           IncI_unscaled = round(IncI*popsize/10^4)) %>% 
    rename(dept_id = id, day = time) %>%
    left_join(data_SEIR_monolix, by = c("dept_id", "day"))
  
  # EpiEstim only
  res_EE <- EpiEstim_only_fun(data_for_est = data_x, 
                              Inc_name = "IncI_unscaled", 
                              meansi = 10.1, stdsi = 8.75) %>%
    mutate(model = sim_res_slow_names[j]) %>%
    left_join(data_x, by = c("dept_id", "day"))
  
  res_EE_slow_NPI_list3[[j]] <- res_EE
  
  # regression
  reg_res <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = res_EE)
  
  # coefficients
  coefs_Rt_reg <- coefficients(reg_res)$dept_id
  
  confint_Rt_reg <- data.frame(confint(reg_res, method="Wald"))[-c(1:2), ]
  names(confint_Rt_reg) <- c("CI_LL", "CI_UL")
  
  
  coefs_Rt_reg <- coefs_Rt_reg %>%
    dplyr::select(-1) %>%
    unique() %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
    cbind(confint_Rt_reg[-1,]) %>%
    mutate(parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("Lockdown 1", "Barrier gestures")), 
           model = sim_res_slow_names[j])
  
  res_reg_slow_NPI_list3[[j]] <- coefs_Rt_reg
  
  
  # regression fits
  fitted_vals <- exp(fitted(reg_res))
  fitted_df <- res_EE %>%
    filter(!is.na(Rt)) %>%
    mutate(Rt_fitted = fitted_vals) %>%
    select(dept_id, day, Rt_SEIRAHD, Rt, CI_LL, CI_UL, Rt_fitted, model, lockdown1, BG1)
  
  fits_reg_slow_NPI_list3[[j]] <- fitted_df
}

stopCluster(cl)


# EpiEstim plots
res_EE_slow_NPI_df3 <- do.call("rbind.data.frame", res_EE_slow_NPI_list3)  %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94))) %>%
  mutate(model = str_replace_all(model, "_", " "))

ggplot(res_EE_slow_NPI_df3 %>% filter(dept_id %in% selected_depts),
       aes(x = day, y = Rt, col = model)) + 
  geom_line() +
  geom_ribbon(aes(ymax = CI_UL, ymin = CI_LL, fill = model), alpha = 0.4, color = NA) +
  geom_line(aes(y = Rt_SEIRAHD), linetype = "dashed") + 
  geom_xsideline(aes(y = IncI_unscaled)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_wrap(~dept_id2) +
  scale_x_continuous(breaks = seq(0, 120, 20)) + 
  labs(title = "Comparison Rt estimated by EpiEstim and real underlying Rt", 
       col = "", fill = "") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + 
  theme(legend.position = "bottom")  + 
  theme(legend.text = element_text(family = "serif", size = 12, hjust = 0), 
        legend.title = element_text(family = "serif", size = 13),
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))


# regression plots
res_reg_slow_NPI_df3 <- do.call("rbind.data.frame", res_reg_slow_NPI_list3) %>%
  mutate(true_value = ifelse(parameter == "Lockdown 1", -1.45, -0.8), 
         bias = abs(true_value - value), 
         rel_bias = abs(true_value - value)/abs(true_value)*100) %>%
  mutate(model = str_replace_all(model, "_", " "))

ggplot(res_reg_slow_NPI_df3, aes(ymin = CI_LL, ymax = CI_UL, x = model, 
                                y = value, col = model)) + 
  geom_pointrange(position = position_dodge(width = 1)) +  
  facet_wrap(~parameter, ncol = 1, scale = "free_y") +
  geom_hline(aes(yintercept = true_value), linetype = "dashed", col = "black", linewidth = 0.8) + 
  labs(title = "Comparison 2 step regression: slower NPI onset", col = "",
       x = "", y = "coefficient value") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.text = element_text(family = "serif", size = 12, hjust = 0), 
        legend.title = element_text(family = "serif", size = 13),
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13))


# regression fits
fits_reg_slow_NPI_df3 <- do.call("rbind.data.frame", fits_reg_slow_NPI_list3) 
fits_reg_plot_slow3 <- res_EE_slow_NPI_df3 %>%
  select(-c(Rt, Rt_SEIRAHD, CI_LL, CI_UL)) %>%
  left_join(., fits_reg_slow_NPI_df3 %>% 
            mutate(model = str_replace_all(model, "_", " ")), 
                   by = c("dept_id", "day", "model")) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94)))


ggplot(fits_reg_plot_slow3 %>% filter(dept_id %in% selected_depts),
       aes(x = day, y = Rt, col = model)) + 
  geom_line(aes(linetype = "EpiEstim Rt"), linewidth = 0.7) +
  geom_ribbon(aes(ymax = CI_UL, ymin = CI_LL, fill = model), alpha = 0.4, color = NA) +
  geom_line(aes(y = Rt_SEIRAHD, linetype = "True Rt"), linewidth = 0.7) + 
  geom_line(aes(y = Rt_fitted, linetype = "Regression fit Rt"), linewidth = 0.7) + 
  geom_xsideline(aes(y = IncI_unscaled)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_wrap(~dept_id2) +
  scale_x_continuous(breaks = seq(0, 120, 20)) + 
  scale_linetype_manual(values = c("True Rt" = "solid", 
                                   "EpiEstim Rt" = "dotted", 
                                   "Regression fit Rt" = "dashed")) + 
  labs(#title = expression(R[t] ~ "estimated by EpiEstim and real underlying" ~ R[t]), 
       col = "", fill = "", linetype = "", x = "Day") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + 
  theme(legend.position = "bottom", 
        legend.text = element_text(family = "serif", size = 12, hjust = 0), 
        legend.title = element_text(family = "serif", size = 13),
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))


ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/fits_slow_NPIs.jpeg", 
       dpi = 400, width = 14, height = 8.5)


# metrics table
metrics_table_slow3 <- res_reg_slow_NPI_df3 %>%
  mutate(reg_coef = paste0(round(value, 3), " [", 
                           round(CI_LL, 3), ", ", round(CI_UL, 3), "]"), 
         bias = as.character(round(bias, 3)), 
         rel_bias = as.character(round(rel_bias, 2))) %>%
  select(-c(value, CI_LL, CI_UL, true_value)) %>%
  pivot_longer(cols = c(reg_coef, bias, rel_bias), 
               names_to = "metric", 
               values_to = "value") %>%
  pivot_wider(names_from = model, values_from = value, 
              names_sort = TRUE) %>%
  mutate(metric = factor(metric, levels = c("reg_coef", "bias", "rel_bias"))) %>%
  arrange (parameter, metric) %>%
  mutate(metric = case_when(metric == "bias" ~ "absolute bias", 
                            metric == "reg_coef" ~ "regression coefficient", 
                            metric == "rel_bias" ~ "relative bias (%)"))


metrics_table_slow3 %>% 
  ungroup() %>%
  dplyr::select(-parameter) %>%
  kable(digits = 2, format = "html" #, table.attr = "style='width:40%;'"
  ) %>%
  kable_styling(bootstrap_options = "striped") %>%
  pack_rows("NPI 1", 1, 3) %>%
  pack_rows("NPI 2", 4, 6) %>%
  column_spec(1, width_min = "3.8cm")


