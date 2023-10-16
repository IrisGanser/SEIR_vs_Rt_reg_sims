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


dir2 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"
dir3 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once4"

popsize_df <- read.csv(paste0(dir2, "/popsize_df.csv")) %>%
  mutate(dept_id = ifelse(dept_id > 20, dept_id-1, dept_id))

load(paste0(dir2, "/ind_param_list.RData"))


#### load data and calculate Rt ####
data_list_Simulx <- loadRData(paste(dir2, "sim_res_Simulx_all_list2.RData", sep = "/"))

dataset1_Simulx <- data_list_Simulx[[1]] %>%
  left_join(popsize_df, by = c("id" = "dept_id")) %>%
  left_join(ind_param_list[[1]], by = "id") %>%
  mutate(Rt_SEIRAHD = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
         IncI_unscaled = round(IncI*popsize/10^4), 
         IncH_unscaled = round(IncH*popsize/10^4), 
         lockdown1 = ifelse(between(time, 16, 70), 1, 0), 
                BG1 = ifelse(time > 70, 1, 0)) %>% 
  rename(dept_id = id, day = time)

ggplot(dataset1_Simulx, aes(x = day, y = Rt_SEIRAHD, group = dept_id)) + 
  geom_line() +
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  labs(title = "Rt from SEIRAHD model", y = "Rt") + 
  theme_bw()


data_ABM_rm_cov <- read.csv(paste0(dir3, "/data_covasim_rm5_Rt_1.csv")) %>%
  filter(day > 15) %>%
  mutate(day = day - 15)
data_ABM_hybrid_cov <- read.csv(paste0(dir3, "/data_covasim_hybrid5_Rt_1.csv")) %>%
  filter(day > 15) %>%
  mutate(day = day - 15)


ggplot(data_ABM_rm_cov, aes(x = day, y = Rt, group = dept_id)) + 
  geom_line() +
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  labs(title = "Rt from ABM model", y = "Rt") + 
  theme_bw()


#### 2 step reg ####
cl <- makeCluster(6)
registerDoParallel(cl)

res_reg_Simulx_I <- EpiEstim_reg_fun(data_for_est = dataset1_Simulx, 
                                   Inc_name = "IncI_unscaled", 
                                   rep_num = 1, 
                                   meansi = 10.1, stdsi = 8.75)

res_reg_Simulx_H <- EpiEstim_reg_fun(data_for_est = dataset1_Simulx, 
                                     Inc_name = "IncI_unscaled", 
                                     rep_num = 1, 
                                     lag_NPIs = TRUE, lag_days = 5, 
                                     meansi = 10.1, stdsi = 8.75)

res_reg_ABM_rm <- EpiEstim_reg_fun(data_for_est = data_ABM_rm_cov, 
                                   Inc_name = "IncI", 
                                   rep_num = 1, 
                                   meansi = 7.68, stdsi = 4.17)

res_reg_ABM_hybrid <- EpiEstim_reg_fun(data_for_est = data_ABM_hybrid_cov, 
                                  Inc_name = "IncI", 
                                  rep_num = 1, 
                                  meansi = 7.68, stdsi = 4.17)

stopCluster(cl)

reg_res_df <- bind_rows(res_reg_Simulx_I %>% mutate(model = "Simulx IncI"), 
                        res_reg_Simulx_H %>% mutate(model = "Simulx IncH"), 
                        res_reg_ABM_rm %>% mutate(model = "ABM rm IncI"),
                        res_reg_ABM_hybrid %>% mutate(model = "ABM hybrid IncI")) %>%
  mutate(parameter = factor(parameter, 
                            levels = c("Lockdown 1", "Barrier gestures"),
                            labels = c("NPI 1", "NPI 2")),
         true_value = ifelse(parameter == "NPI 1", -1.45, -0.5))


br_palette <- diverging_hcl("Blue-Red", n = 20)
pg_palette <- diverging_hcl("Purple-Green", n = 20)
plot_cols <- c(pg_palette[c(18, 14)], br_palette[c(1, 5)])

ggplot(reg_res_df, aes(x = 1, y = value, ymin = CI_LL, ymax = CI_UL, col = model)) + 
  geom_pointrange(position = position_dodge(width = 1))  +
  geom_hline(data = data.frame(yval = -1.45, parameter = "NPI 1"), aes(yintercept = yval),
             linetype = "dashed", col = "black", linewidth = 0.8) +
  geom_hline(data = data.frame(yval = -0.5, parameter = "NPI 2"), aes(yintercept = yval),
             linetype = "dashed", col = "black", linewidth = 0.8) + 
  facet_wrap(~parameter) +
  labs(x = "", y = "coefficient value", 
       title = "Regression estimates") + 
  scale_color_manual(values = plot_cols) +
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) 



#### only EpiEstim ####
cl <- makeCluster(6)
registerDoParallel(cl)

res_EpiEstim_Simulx_I <- EpiEstim_only_fun(data_for_est = dataset1_Simulx, 
                                     Inc_name = "IncI_unscaled", 
                                     meansi = 10.1, stdsi = 8.75)

res_EpiEstim_Simulx_H <- EpiEstim_only_fun(data_for_est = dataset1_Simulx, 
                                     Inc_name = "IncI_unscaled", 
                                     meansi = 10.1, stdsi = 8.75)

res_EpiEstim_ABM_rm <- EpiEstim_only_fun(data_for_est = data_ABM_rm_cov, 
                                  Inc_name = "IncI", 
                                  meansi = 7.41, stdsi = 3.96)

res_EpiEstim_ABM_hybrid <- EpiEstim_only_fun(data_for_est = data_ABM_hybrid_cov, 
                                    Inc_name = "IncI", 
                                    meansi =  7.31, stdsi = 3.81)

stopCluster(cl)

Rt_comp_df <- bind_rows(res_EpiEstim_Simulx_I %>%
                          left_join(dataset1_Simulx %>% rename(Rt_real = Rt_SEIRAHD), 
                                    by = c("dept_id", "day")) %>%
                          mutate(model = "Simulx IncI"), 
                        res_EpiEstim_Simulx_H %>%
                          left_join(dataset1_Simulx %>% rename(Rt_real = Rt_SEIRAHD), 
                                    by = c("dept_id", "day")) %>%
                          mutate(model = "Simulx IncH"),
                        res_EpiEstim_ABM_rm %>%
                          left_join(data_ABM_rm_cov %>% rename(Rt_real = Rt), 
                                    by = c("dept_id", "day")) %>%
                          mutate(model = "ABM rm IncI"),
                        res_EpiEstim_ABM_hybrid %>%
                          left_join(data_ABM_hybrid_cov %>% rename(Rt_real = Rt), 
                                    by = c("dept_id", "day")) %>%
                          mutate(model = "ABM hybrid IncI"))

ggplot(Rt_comp_df %>% filter(dept_id %in% c(1, 4, 9, 13)),
       aes(x = day, y = Rt)) + 
  geom_line(aes(col = "EpiEstim")) +
  geom_ribbon(aes(ymax = CI_UL, ymin = CI_LL, fill = "EpiEstim"), alpha = 0.5, show.legend = FALSE) +
  geom_line(aes(y = Rt_real, col = "Real Rt")) + 
  facet_grid(rows = vars(model), cols = vars(dept_id)) + 
  labs(title = expression("Comparison"~  R[t]~ "estimated by EpiEstim and real underlying" ~ R[t]), 
       col = "", fill = "", y = expression(R[t])) +
  theme_bw() +
  scale_color_brewer(palette = "Set1", 
                     labels = c(expression("EpiEstim" ~R[t]), expression("Real" ~R[t]))) +
  scale_fill_brewer(palette = "Set1", 
                    labels = c(expression("EpiEstim" ~R[t]), expression("Real" ~R[t]))) + 
  theme(legend.text = element_text(hjust = 0))


save(Rt_comp_df, reg_res_df, file = "reg_Rt_example.RData")


#### regressions with correct Rt ####
reg_Rt_known_Simulx <- lmer(log(Rt_SEIRAHD) ~ lockdown1 + BG1 + (1|dept_id), data = dataset1_Simulx)
reg_Rt_known_ABM_rm <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_ABM_rm_cov)
reg_Rt_known_ABM_hybrid <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_ABM_hybrid_cov)

reg_models_list <- list(reg_Rt_known_Simulx, reg_Rt_known_ABM_rm, reg_Rt_known_ABM_hybrid)
model_names <- c("Simulx", "ABM rm", "ABM hybrid")

for(i in 1:length(reg_models_list)){
  residuals <- residuals(reg_models_list[[i]])
  fitted_vals <- fitted(reg_models_list[[i]])
  
  plot(exp(fitted_vals), residuals, 
       xlab = "Fitted Rt", ylab = "Residual", 
       main = paste("Diagnostics True Rt regression", model_names[i]))
}

reg_Rt_known_comp_list <- list()
for(i in 1:length(reg_models_list)){
  
  coefs_Rt_reg <- coefficients(reg_models_list[[i]])$dept_id
  
  confint_Rt_reg <- data.frame(confint(reg_models_list[[i]], method="Wald"))[-c(1:2), ]
  names(confint_Rt_reg) <- c("CI_LL", "CI_UL")
  
  
  reg_Rt_known_comp_list[[i]] <- coefs_Rt_reg %>%
    dplyr::select(-1) %>%
    unique() %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
    cbind(confint_Rt_reg[-1,]) %>%
    mutate(parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("Lockdown 1", "Barrier gestures"))) %>%
    mutate(model = model_names[i])
}

reg_Rt_known_comp <- do.call("rbind.data.frame", reg_Rt_known_comp_list) %>%
  mutate(true_value = ifelse(parameter == "Lockdown 1", -1.45, -0.5))

ggplot(reg_Rt_known_comp, aes(x = 1, y = value, ymin = CI_LL, ymax = CI_UL, col = model)) + 
  geom_pointrange(position = position_dodge(width = 1)) + 
  facet_wrap(~parameter, scales = "free_x") +
  geom_line(aes(y = true_value, x = seq(0.2, 2, length.out = 6)), 
            linetype = "dashed", col = "black", linewidth = 0.8) + 
  labs(x = "", y = "coefficient value", 
       title = "Regression estimates with correct Rt") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  scale_color_brewer(palette = "Dark2")


# bias table
metric_table_reg_Rt_known <- reg_Rt_known_comp %>%
  mutate(parameter = ifelse(parameter == "Lockdown 1", "NPI 1", "NPI 2"), 
         coverage = ifelse(between(true_value, CI_LL, CI_UL), "yes", "no"), 
         bias = as.character(round(true_value - value, 2)), 
         relbias = as.character(round(abs(true_value - value)/abs(true_value)*100, 2))) %>%
  select(parameter, model, coverage, bias, relbias) %>%
  unique() %>%
  pivot_longer(cols = c(coverage, bias, relbias), 
               names_to = "metric", 
               values_to = "value") %>%
  pivot_wider(names_from = model, values_from = value) %>%
  mutate(metric = factor(metric, levels = c("coverage", "bias", "relbias"))) %>%
  arrange (parameter, metric) %>%
  mutate(metric = case_when(metric == "bias" ~ "absolute bias", 
                            metric == "coverage" ~ "CI coverage", 
                            metric == "relbias" ~ "relative bias (%)"))

metric_table_reg_Rt_known %>% 
  select(-parameter) %>%
  kable(digits = 2, format = "html") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
  pack_rows("NPI 1", 1, 3) %>%
  pack_rows("NPI 2", 4, 6)



# regression fits
data_list <- list(dataset1_Simulx %>% rename(Rt = Rt_SEIRAHD), data_ABM_rm_cov, data_ABM_hybrid_cov)

for(i in 1:length(reg_models_list)){
  data <- data.frame(dept_id = data_list[[i]]$dept_id[data_list[[i]]$Rt > 0],
                     day = data_list[[i]]$day[data_list[[i]]$Rt > 0],
                     true_Rt = data_list[[i]]$Rt[data_list[[i]]$Rt > 0], 
                     fitted_Rt = exp(fitted(reg_models_list[[i]])), 
                     model = model_names[i])
  
  plot <- ggplot(data %>% filter(dept_id %in% 1:12), aes(x = day, y = true_Rt)) + 
    geom_line(aes(col = "True Rt")) + 
    geom_line(aes(y = fitted_Rt, col = "Fitted Rt")) + 
    facet_wrap(~dept_id) + 
    scale_color_brewer(palette = "Dark2") + 
    theme_bw() + 
    labs(y = "Rt", col = "",
         title = paste("Rt comparison", model_names[i]))
  
  print(plot)
  
  assign(paste0("reg_fits_", model_names[i]), data)
}

reg_fits_all <- bind_rows(reg_fits_Simulx, `reg_fits_ABM rm`, `reg_fits_ABM hybrid`)

ggplot(reg_fits_all %>% filter(dept_id %in% c(1, 4, 9, 13)), 
       aes(x = day, y = true_Rt)) + 
  geom_line(aes(col = "Real Rt")) + 
  geom_line(aes(y = fitted_Rt, col = "Fitted Rt")) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(rows = vars(model), cols = vars(dept_id)) +
  scale_color_brewer(palette = "Set1", 
                     labels = c(expression("Fitted" ~R[t]), expression("Real" ~R[t]))) + 
  theme_bw() + 
  labs(y = expression(R[t]), col = "",
       title = expression("Fits from regression with known"~ R[t])) +
  theme(legend.text = element_text(hjust = 0))


save(reg_Rt_known_comp, metric_table_reg_Rt_known, reg_fits_all,
     file = "res_reg_Rt_known.RData")

