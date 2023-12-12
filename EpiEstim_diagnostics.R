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
library(ggridges)
library(kableExtra)

setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/Rt_trajectories")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/deSolve_SIR_model_function.R")
source("~/PhD/COVID_France/Dropbox_iris_covid/departement/Données_SPF/Data/data_functions.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")


dir2 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"
dir3 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once7"
dir4 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once8"



# load data and calculate Rt ----------------------------------------------
popsize_df <- read.csv(paste0(dir2, "/popsize_df.csv")) %>%
  mutate(dept_id = ifelse(dept_id > 20, dept_id-1, dept_id))

indparams_Simulx2_1 <- read.table(paste0(dir2, "/ind_params/ind_2params_1.txt"), 
                                  header = TRUE, sep = " ")
indparams_Simulx4_1 <- read.table(paste0(dir2, "/ind_params/ind_2params_new4_1.txt"), 
                                  header = TRUE, sep = " ")


data_list_Simulx_new2 <- loadRData(paste(dir2, "sim_res_Simulx_2params_new2_list.RData", sep = "/"))
data_list_Simulx_new4 <- loadRData(paste(dir2, "sim_res_Simulx_2params_new4_list.RData", sep = "/"))

dataset1_Simulx_new2 <- data_list_Simulx_new2[[1]] %>%
  left_join(popsize_df, by = c("id" = "dept_id")) %>%
  left_join(indparams_Simulx2_1, by = "id") %>%
  mutate(Rt = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
         IncI_unscaled = IncI*popsize/10^4, 
         IncH_unscaled = IncH*popsize/10^4,
         lockdown1 = ifelse(between(time, 16, 70), 1, 0), 
         BG1 = ifelse(time > 70, 1, 0)) %>% 
  rename(dept_id = id, day = time)

true_Rt_df_Simulx_new2 <- dataset1_Simulx_new2 %>% dplyr::select(dept_id, day, Rt) %>% 
  rename(Rt_real = Rt)


dataset1_Simulx_new4 <- data_list_Simulx_new4[[1]] %>%
  left_join(popsize_df, by = c("id" = "dept_id")) %>%
  left_join(indparams_Simulx4_1, by = "id") %>%
  mutate(Rt = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
         IncI_unscaled = IncI*popsize/10^4, 
         IncH_unscaled = IncH*popsize/10^4,
         lockdown1 = ifelse(between(time, 45, 99), 1, 0), 
         BG1 = ifelse(time > 99, 1, 0)) %>% 
  rename(dept_id = id, day = time)

true_Rt_df_Simulx_new4 <- dataset1_Simulx_new4 %>% dplyr::select(dept_id, day, Rt) %>% 
  rename(Rt_real = Rt)



data_ABM_rm_cov <- read.csv(paste0(dir3, "/data_covasim_rm7_Rt_1.csv")) %>%
  filter(day > 15) %>%
  mutate(day = day - 15)
data_ABM_hybrid_cov <- read.csv(paste0(dir3, "/data_covasim_hybrid7_Rt_1.csv")) %>%
  filter(day > 15) %>%
  mutate(day = day - 15)


data_ABM_rm_cov_long <- read.csv(paste0(dir3, "/data_covasim_rm7_Rt_1.csv")) 
data_ABM_hybrid_cov_long <- read.csv(paste0(dir3, "/data_covasim_hybrid7_Rt_1.csv"))

true_Rt_df_ABM_rm_long <- data_ABM_rm_cov_long %>% 
  dplyr::select(dept_id, day, Rt) %>% 
  rename(Rt_real = Rt)
true_Rt_df_ABM_hybrid_long <- data_ABM_hybrid_cov_long %>% 
  dplyr::select(dept_id, day, Rt) %>% 
  rename(Rt_real = Rt)


set.seed(123)
selected_depts <- sample(1:94, size = 16, replace = FALSE)


# known Rt regressions ------------------------------------------------

# run regressions
reg_Rt_known_Simulx2 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = dataset1_Simulx_new2)
reg_Rt_known_Simulx4 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = dataset1_Simulx_new4)
reg_Rt_known_ABM_rm <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_ABM_rm_cov)
reg_Rt_known_ABM_hybrid <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_ABM_hybrid_cov)
reg_Rt_known_ABM_rm_long <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_ABM_rm_cov_long)
reg_Rt_known_ABM_hybrid_long <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_ABM_hybrid_cov_long)

reg_models_list <- list(reg_Rt_known_Simulx2, reg_Rt_known_Simulx4,
                        reg_Rt_known_ABM_rm, reg_Rt_known_ABM_hybrid,
                        reg_Rt_known_ABM_rm_long, reg_Rt_known_ABM_hybrid_long)
model_names <- c("SEIRAHD", "SEIRAHD long", "random mixing ABM", "hybrid ABM", 
                 "random mixing ABM long", "hybrid ABM long")


# initialize result list
reg_Rt_known_comp_list <- list()

# summarize results
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
                              labels = c("NPI 1", "NPI 2"))) %>%
    mutate(model = model_names[i])
}

reg_Rt_known_comp <- do.call("rbind.data.frame", reg_Rt_known_comp_list) %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.5))

ggplot(reg_Rt_known_comp, aes(x = model, y = value, ymin = CI_LL, ymax = CI_UL, col = model)) + 
  geom_pointrange(position = position_dodge(width = 1)) + 
  facet_wrap(~parameter, scales = "free_x") +
  geom_hline(data = data.frame(true_value = c(-1.45, -0.5), parameter = c("NPI 1", "NPI 2")), 
             aes(yintercept = true_value), linetype = "dashed") + 
  labs(x = "", y = "coefficient value", 
       title = "Regression estimates with correct Rt") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  scale_color_brewer(palette = "Dark2")


# fits
data_list <- list(dataset1_Simulx_new2, dataset1_Simulx_new4, data_ABM_rm_cov, 
                  data_ABM_hybrid_cov, data_ABM_rm_cov_long, data_ABM_hybrid_cov_long)
fit_data_list <- list()

for(i in 1:length(reg_models_list)){
  
  fitted <- exp(fitted(reg_models_list[[i]]))
  
  fit_data_list[[i]] <- data_list[[i]] %>%
    mutate(fitted_Rt = fitted, 
           model = model_names[i])
}

reg_fits_all <- do.call("bind_rows", fit_data_list) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)))


ggplot(reg_fits_all %>% filter(dept_id %in% c(selected_depts[1:4])), 
       aes(x = day, y = Rt)) + 
  geom_line(aes(col = "Real Rt")) + 
  geom_line(aes(y = fitted_Rt, col = "Fitted Rt")) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(rows = vars(model), cols = vars(dept_id2)) +
  scale_color_brewer(palette = "Set1", 
                     labels = c(expression("Fitted" ~R[t]), expression("Real" ~R[t]))) + 
  theme_bw() + 
  labs(y = expression(R[t]), col = "",
       title = expression("Fits from regression with known"~ R[t])) +
  theme(legend.text = element_text(hjust = 0))


# bias table
metric_table_reg_Rt_known <- reg_Rt_known_comp %>%
  mutate(coverage = ifelse(between(true_value, CI_LL, CI_UL), "yes", "no"), 
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




# Fits EpiEstim + reg -----------------------------------------------------
start_days2 <- 2:119
end_days2 <- 4:121

start_days4 <- 2:148
end_days4 <- 4:150

cl <- makeCluster(10)
registerDoParallel(cl)

Rt_comp_res2 <- Rt_calc_fun(Inc_df = dataset1_Simulx_new2, id_col_name = "dept_id", time_col_name = "day", 
                            Inc_col_name = "IncI_unscaled", model_name = "lag 0", 
                            Rt_ref_df = true_Rt_df_Simulx_new2, tstart = start_days2, tend = end_days2,
                            meansi = 10.1, stdsi = 8.75, Rt_prior = 1)

Rt_comp_res4 <- Rt_calc_fun(Inc_df = dataset1_Simulx_new4, id_col_name = "dept_id", time_col_name = "day", 
                            Inc_col_name = "IncI_unscaled", model_name = "lag 0", 
                            Rt_ref_df = true_Rt_df_Simulx_new4, tstart = start_days4, tend = end_days4, 
                            meansi = 10.1, stdsi = 8.75, Rt_prior = 1)
stopCluster(cl)

# summary SEIRAHD 2 dataset
reg_data2 <- Rt_comp_res2$Rt_comp %>%
  mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
         BG1 = ifelse(day > 70, 1, 0)) %>%
  filter(!is.na(Rt))

reg_res2 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data2)

fitted_vals2 <- exp(fitted(reg_res2))
fitted_df2 <- reg_data2 %>%
  filter(!is.na(Rt)) %>%
  mutate(Rt_fitted = fitted_vals2)

# summary SEIRAHD 4 dataset
reg_data4 <- Rt_comp_res4$Rt_comp %>%
  mutate(lockdown1 = ifelse(between(day, 45, 99), 1, 0), 
         BG1 = ifelse(day > 99, 1, 0)) %>%
  filter(!is.na(Rt))

reg_res4 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data4)

fitted_vals4 <- exp(fitted(reg_res4))
fitted_df4 <- reg_data4 %>%
  filter(!is.na(Rt)) %>%
  mutate(Rt_fitted = fitted_vals4)


# summary SEIRAHD 4 dataset with first days cut off
reg_data_m10 <- Rt_comp_res4$Rt_comp %>%
  mutate(lockdown1 = ifelse(between(day, 45, 99), 1, 0), 
         BG1 = ifelse(day > 99, 1, 0)) %>%
  filter(day > 13)

reg_res_m10 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data_m10)

fitted_vals_m10 <- exp(fitted(reg_res_m10))
fitted_df_m10 <- reg_data_m10 %>%
  filter(day > 13) %>%
  mutate(Rt_fitted = fitted_vals_m10)


rect_cols <- sequential_hcl(5, palette = "BluYl")

ggplot(fitted_df2 %>% filter(dept_id %in% selected_depts), aes(x = day, y = Rt_fitted)) + 
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = 43, y = Inf, label = "NPI 1", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  annotate("label", x = 96, y = Inf, label =  "NPI 2", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  geom_line(aes(linetype = "Regression fit Rt")) +
  geom_line(aes(y = Rt_real, linetype = "Real Rt")) +
  geom_line(aes(y = Rt, linetype = "EpiEstim Rt")) +
  labs(linetype = "") + 
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 121)) + 
  labs(y = expression(R[t]), title = "Simulx infections, short NPI-free period, no days cut off, no NPI lag") +
  facet_wrap(~dept_id) + 
  theme_bw()

ggplot(fitted_df4 %>% filter(dept_id %in% selected_depts), aes(x = day, y = Rt_fitted)) + 
  annotate("rect", xmin = 45, xmax = 100, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 100, xmax = 151, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = 75, y = Inf, label = "NPI 1", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  annotate("label", x = 125, y = Inf, label =  "NPI 2", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  geom_line(aes(linetype = "Regression fit Rt")) +
  geom_line(aes(y = Rt_real, linetype = "Real Rt")) +
  geom_line(aes(y = Rt, linetype = "EpiEstim Rt")) +
  labs(linetype = "") + 
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 151)) + 
  labs(y = expression(R[t]), title = "Simulx infections, long NPI-free period, no days cut off, no NPI lag") +
  facet_wrap(~dept_id) + 
  theme_bw()


ggplot(fitted_df_m10 %>% filter(dept_id %in% selected_depts), aes(x = day, y = Rt_fitted)) + 
  annotate("rect", xmin = 45, xmax = 100, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 100, xmax = 151, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = 75, y = Inf, label = "NPI 1", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  annotate("label", x = 125, y = Inf, label =  "NPI 2", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  geom_line(aes(linetype = "Regression fit Rt")) +
  geom_line(aes(y = Rt_real, linetype = "Real Rt")) +
  geom_line(aes(y = Rt, linetype = "EpiEstim Rt")) +
  labs(linetype = "") + 
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 151)) + 
  labs(y = expression(R[t]), title = "Simulx infections, long NPI-free period, first 10 days cut off, no NPI lag") +
  facet_wrap(~dept_id) + 
  theme_bw()


# bias assessment
bias_first_days2 <- fitted_df2 %>%
  group_by(dept_id) %>%
  filter(between(day, 8, 15)) %>%
  summarize(bias = mean(Rt_fitted - Rt_real))

bias_first_days4 <- fitted_df4 %>%
  group_by(dept_id) %>%
  filter(between(day, 8, 44)) %>%
  summarize(bias = mean(Rt_fitted - Rt_real))

bias_first_days_m10 <- fitted_df_m10 %>%
  group_by(dept_id) %>%
  filter(between(day, 8, 44)) %>%
  summarize(bias = mean(Rt_fitted - Rt_real))

bias_first_days <- bias_first_days2 %>% mutate(model = "short") %>%
  bind_rows(bias_first_days4 %>% mutate(model = "long")) %>%
  bind_rows(bias_first_days_m10 %>% mutate(model = "long m10")) %>%
  mutate(model = as.factor(model))

ggplot(bias_first_days, aes(x = as.factor(dept_id), y = bias, fill = model)) + 
  geom_col(position = position_dodge(width = 0.8)) +
  scale_fill_brewer(palette = "Dark2")
