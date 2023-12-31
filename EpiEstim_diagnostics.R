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
library(ggside)

setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/Rt_trajectories")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/deSolve_SIR_model_function.R")
source("~/PhD/COVID_France/Dropbox_iris_covid/departement/Données_SPF/Data/data_functions.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")


dir2 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"
dir3 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once7"
dir4 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once8"
dir10 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once10"
dir11 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once11"


# load data and calculate Rt ----------------------------------------------
popsize_df <- read.csv(paste0(dir2, "/popsize_df.csv")) %>%
  mutate(dept_id = ifelse(dept_id > 20, dept_id-1, dept_id))

indparams_Simulx2_1 <- read.table(paste0(dir2, "/ind_params/ind_2params_1.txt"), 
                                  header = TRUE, sep = " ")
indparams_Simulx3_1 <- read.table(paste0(dir2, "/ind_params/ind_2params_new3_1.txt"), 
                                  header = TRUE, sep = " ")
indparams_Simulx4_1 <- read.table(paste0(dir2, "/ind_params/ind_2params_new4_1.txt"), 
                                  header = TRUE, sep = " ")
indparams_Simulx4_1 <- read.table(paste0(dir2, "/ind_params/ind_2params_new5_1.txt"), 
                                  header = TRUE, sep = " ")


data_list_Simulx_new2 <- loadRData(paste(dir2, "sim_res_Simulx_2params_new2_list.RData", sep = "/"))
data_list_Simulx_new3 <- loadRData(paste(dir2, "sim_res_Simulx_2params_new3_list.RData", sep = "/"))
data_list_Simulx_new4 <- loadRData(paste(dir2, "sim_res_Simulx_2params_new4_list.RData", sep = "/"))
data_list_Simulx_new5 <- loadRData(paste(dir2, "sim_res_Simulx_2params_new5_list.RData", sep = "/"))


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


dataset1_Simulx_new3 <- data_list_Simulx_new3[[1]] %>%
  left_join(popsize_df, by = c("id" = "dept_id")) %>%
  left_join(indparams_Simulx3_1, by = "id") %>%
  mutate(Rt = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
         IncI_unscaled = IncI*popsize/10^4, 
         IncH_unscaled = IncH*popsize/10^4,
         lockdown1 = ifelse(between(time, 16, 70), 1, 0), 
         BG1 = ifelse(time > 70, 1, 0)) %>% 
  rename(dept_id = id, day = time)

true_Rt_df_Simulx_new3 <- dataset1_Simulx_new3 %>% dplyr::select(dept_id, day, Rt) %>% 
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


dataset1_Simulx_new5 <- data_list_Simulx_new5[[1]] %>%
  left_join(popsize_df, by = c("id" = "dept_id")) %>%
  left_join(indparams_Simulx5_1, by = "id") %>%
  mutate(Rt = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
         IncI_unscaled = IncI*popsize/10^4, 
         IncH_unscaled = IncH*popsize/10^4, 
         lockdown1 = ifelse(between(time, 45, 99), 1, 0), 
         BG1 = ifelse(time > 99, 1, 0)) %>% 
  rename(dept_id = id, day = time)

true_Rt_df_Simulx_new5 <- dataset1_Simulx_new5 %>% dplyr::select(dept_id, day, Rt) %>% 
  rename(Rt_real = Rt)



data_ABM_rm_cov <- read.csv(paste0(dir3, "/data_covasim_rm7_Rt_1.csv")) %>%
  filter(day > 15) %>%
  mutate(day = day - 15)
data_ABM_hybrid_cov <- read.csv(paste0(dir3, "/data_covasim_hybrid7_Rt_1.csv")) %>%
  filter(day > 15) %>%
  mutate(day = day - 15)


data_ABM_rm_cov10 <- read.csv(paste0(dir10, "/data_covasim_rm10_Rt_1.csv")) %>%
  filter(day > 29) %>%
  mutate(day = day - 29)
data_ABM_hybrid_cov10 <- read.csv(paste0(dir10, "/data_covasim_hybrid10_Rt_1.csv")) %>%
  filter(day > 29) %>%
  mutate(day = day - 29)

data_ABM_rm_cov_long <- read.csv(paste0(dir3, "/data_covasim_rm7_Rt_1.csv")) 
data_ABM_hybrid_cov_long <- read.csv(paste0(dir3, "/data_covasim_hybrid7_Rt_1.csv"))

data_ABM_rm_cov10_long <- read.csv(paste0(dir10, "/data_covasim_rm10_Rt_1.csv")) 
data_ABM_hybrid_cov10_long <- read.csv(paste0(dir10, "/data_covasim_hybrid10_Rt_1.csv"))


true_Rt_df_ABM_rm_long <- data_ABM_rm_cov_long %>% 
  dplyr::select(dept_id, day, Rt) %>% 
  rename(Rt_real = Rt)
true_Rt_df_ABM_hybrid_long <- data_ABM_hybrid_cov_long %>% 
  dplyr::select(dept_id, day, Rt) %>% 
  rename(Rt_real = Rt)


true_Rt_df_ABM_rm_long10 <- data_ABM_rm_cov10_long %>% 
  dplyr::select(dept_id, day, Rt) %>% 
  rename(Rt_real = Rt)
true_Rt_df_ABM_hybrid_long10 <- data_ABM_hybrid_cov10_long %>% 
  dplyr::select(dept_id, day, Rt) %>% 
  rename(Rt_real = Rt)

set.seed(123)
selected_depts <- sample(1:94, size = 16, replace = FALSE)


# known Rt regressions ------------------------------------------------

# run regressions
reg_Rt_known_Simulx2 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = dataset1_Simulx_new2)
reg_Rt_known_Simulx3 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = dataset1_Simulx_new3)
reg_Rt_known_Simulx4 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = dataset1_Simulx_new4)
reg_Rt_known_Simulx5 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = dataset1_Simulx_new5)

reg_Rt_known_ABM_rm <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_ABM_rm_cov)
reg_Rt_known_ABM_hybrid <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_ABM_hybrid_cov)
reg_Rt_known_ABM_rm_long <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_ABM_rm_cov_long)
reg_Rt_known_ABM_hybrid_long <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_ABM_hybrid_cov_long)

reg_Rt_known_ABM_rm10 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_ABM_rm_cov10)
reg_Rt_known_ABM_hybrid10 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_ABM_hybrid_cov10)
reg_Rt_known_ABM_rm_long10 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_ABM_rm_cov10_long)
reg_Rt_known_ABM_hybrid_long10 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = data_ABM_hybrid_cov10_long)

reg_models_list <- list(reg_Rt_known_Simulx2, reg_Rt_known_Simulx3, 
                        reg_Rt_known_Simulx4, reg_Rt_known_Simulx5,
                        reg_Rt_known_ABM_rm, reg_Rt_known_ABM_hybrid,
                        reg_Rt_known_ABM_rm_long, reg_Rt_known_ABM_hybrid_long, 
                        reg_Rt_known_ABM_rm10, reg_Rt_known_ABM_hybrid10,
                        reg_Rt_known_ABM_rm_long10, reg_Rt_known_ABM_hybrid_long10)
model_names <- c("SEIRAHD", "SEIRAHD3", "SEIRAHD long", "SEIRAHD5", 
                 "random mixing ABM", "hybrid ABM", 
                 "random mixing ABM long", "hybrid ABM long", 
                 "random mixing ABM 10", "hybrid ABM 10", 
                 "random mixing ABM long 10", "hybrid ABM long 10")


BG1_0.8_models <- c("SEIRAHD3", "random mixing ABM 10", "hybrid ABM 10", 
                    "random mixing ABM long 10", "hybrid ABM long 10")

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
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, 
                             ifelse(model %in% BG1_0.8_models, -0.8, -0.5)))

ggplot(reg_Rt_known_comp %>% filter(model %in% BG1_0.8_models & !grepl("long", model)), 
       aes(x = model, y = value, ymin = CI_LL, ymax = CI_UL, col = model)) + 
  geom_pointrange(position = position_dodge(width = 1)) + 
  facet_wrap(~parameter, scales = "free_x") +
  geom_hline(data = data.frame(true_value = c(-1.45, -0.8), parameter = c("NPI 1", "NPI 2")), 
             aes(yintercept = true_value), linetype = "dashed") + 
  labs(x = "", y = "coefficient value", 
       title = "Regression estimates with correct Rt") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  scale_color_brewer(palette = "Dark2")


# fits
data_list <- list(dataset1_Simulx_new2, dataset1_Simulx_new3, 
                  dataset1_Simulx_new4, dataset1_Simulx_new5, 
                  data_ABM_rm_cov, data_ABM_hybrid_cov,
                  data_ABM_rm_cov_long, data_ABM_hybrid_cov_long, 
                  data_ABM_rm_cov10, data_ABM_hybrid_cov10,
                  data_ABM_rm_cov10_long, data_ABM_hybrid_cov10_long)
fit_data_list <- list()

for(i in 1:length(reg_models_list)){
  
  fitted <- exp(fitted(reg_models_list[[i]]))
  
  fit_data_list[[i]] <- data_list[[i]] %>%
    mutate(fitted_Rt = fitted, 
           model = model_names[i])
}

reg_fits_all <- do.call("bind_rows", fit_data_list) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         IncI_unscaled = ifelse(is.na(IncI_unscaled), IncI, IncI_unscaled))


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


reg_fits_sel <- reg_fits_all %>%
  filter(model %in% BG1_0.8_models & !grepl("long", model)) %>%
  mutate(model = case_when(model == "SEIRAHD3" ~ "SEIRAHD", 
                           grepl("hybrid", model) ~ "hybrid ABM", 
                           grepl("random", model) ~ "random mixing ABM"))

ggplot(reg_fits_sel %>% filter(dept_id %in% selected_depts[c(1:3,5)]), 
       aes(x = day, y = Rt)) + 
  geom_line(aes(col = "Real Rt")) + 
  geom_line(aes(y = fitted_Rt, col = "Fitted Rt")) + 
  geom_xsideline(aes(y = IncI_unscaled)) + 
  ggside(scales = "free_y") + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4)) + 
  facet_grid(rows = vars(model), cols = vars(dept_id2)) +
  scale_color_brewer(palette = "Set1", 
                     labels = c(expression("Fitted" ~R[t]), expression("Real" ~R[t]))) + 
  theme_bw() + 
  labs(y = expression(R[t]), col = "", x = "Day",
       #title = expression("Fits from regression with known"~ R[t])
       ) +
  theme(legend.text = element_text(family = "serif", size = 13, hjust = 0), 
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 11))

ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/reg_fits_Rt_known_3models_wo_title.jpeg", dpi = 400, width = 10, height = 7.5)

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
  filter(model %in% BG1_0.8_models & !grepl("long", model)) %>%
  mutate(model = case_when(model == "SEIRAHD3" ~ "SEIRAHD", 
                           grepl("hybrid", model) ~ "hybrid ABM", 
                           grepl("random", model) ~ "random mixing ABM")) %>%
  pivot_wider(names_from = model, values_from = value) %>%
  mutate(metric = factor(metric, levels = c("coverage", "bias", "relbias"))) %>%
  arrange (parameter, metric) %>%
  mutate(metric = case_when(metric == "bias" ~ "absolute bias", 
                            metric == "coverage" ~ "CI coverage", 
                            metric == "relbias" ~ "relative bias (%)"))

metric_table_reg_Rt_known %>% 
  select(-parameter) %>%
  kable(digits = 2, format = "latex") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
  pack_rows("NPI 1", 1, 3) %>%
  pack_rows("NPI 2", 4, 6)




# Fits EpiEstim + reg -----------------------------------------------------


## Rt estimation -----------------------------------------------------------

cl <- makeCluster(6)
registerDoParallel(cl)

Rt_comp_res2 <- Rt_calc_fun(Inc_df = dataset1_Simulx_new2, id_col_name = "dept_id", time_col_name = "day", 
                            Inc_col_name = "IncI_unscaled", model_name = "lag 0", 
                            Rt_ref_df = true_Rt_df_Simulx_new2,
                            meansi = 10.1, stdsi = 8.75, Rt_prior = 1)

Rt_comp_res3 <- Rt_calc_fun(Inc_df = dataset1_Simulx_new3, id_col_name = "dept_id", time_col_name = "day", 
                            Inc_col_name = "IncI_unscaled", model_name = "lag 0", 
                            Rt_ref_df = true_Rt_df_Simulx_new3,
                            meansi = 10.1, stdsi = 8.75, Rt_prior = 1)

Rt_comp_res4 <- Rt_calc_fun(Inc_df = dataset1_Simulx_new4, id_col_name = "dept_id", time_col_name = "day", 
                            Inc_col_name = "IncI_unscaled", model_name = "lag 0", 
                            Rt_ref_df = true_Rt_df_Simulx_new4,
                            meansi = 10.1, stdsi = 8.75, Rt_prior = 1)


Rt_comp_res5 <- Rt_calc_fun(Inc_df = dataset1_Simulx_new5, id_col_name = "dept_id", time_col_name = "day", 
                            Inc_col_name = "IncI_unscaled", model_name = "lag 0", 
                            Rt_ref_df = true_Rt_df_Simulx_new5,
                            meansi = 10.1, stdsi = 8.75, Rt_prior = 1)
stopCluster(cl)


## Regressions -------------------------------------------------------------
# all with NPIs lagged by 5 days

# summary SEIRAHD 2 dataset
reg_data2 <- Rt_comp_res2$Rt_comp %>%
  mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
         BG1 = ifelse(day > 70, 1, 0)) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         across(lockdown1:BG1, 
                function(x) lag(x, 5, default = 0), 
                .names = "{.col}_lag")) %>%
  filter(!is.na(Rt))

reg_res2 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data2)

fitted_vals2 <- exp(fitted(reg_res2))
fitted_df2 <- reg_data2 %>%
  filter(!is.na(Rt)) %>%
  mutate(Rt_fitted = fitted_vals2)

dataset1_Simulx_new2_fitted <- dataset1_Simulx_new2 %>%
  left_join(fitted_df2 %>% rename(Rt_EE = Rt), by = c("dept_id", "day")) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)),
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94), labels = paste("Region", 1:94)))

# summary SEIRAHD 3 dataset
reg_data3 <- Rt_comp_res3$Rt_comp %>%
  mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
         BG1 = ifelse(day > 70, 1, 0))  %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         across(lockdown1:BG1, 
                function(x) lag(x, 5, default = 0), 
                .names = "{.col}_lag")) %>%
  filter(!is.na(Rt))

reg_res3 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data3)

fitted_vals3 <- exp(fitted(reg_res3))
fitted_df3 <- reg_data3 %>%
  filter(!is.na(Rt)) %>%
  mutate(Rt_fitted = fitted_vals3)

dataset1_Simulx_new3_fitted <- dataset1_Simulx_new3 %>%
  left_join(fitted_df3 %>% rename(Rt_EE = Rt), by = c("dept_id", "day"))%>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94), labels = paste("Region", 1:94)))


# summary SEIRAHD 4 dataset
reg_data4 <- Rt_comp_res4$Rt_comp %>%
  mutate(lockdown1 = ifelse(between(day, 45, 99), 1, 0), 
         BG1 = ifelse(day > 99, 1, 0))  %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         across(lockdown1:BG1, 
                function(x) lag(x, 5, default = 0), 
                .names = "{.col}_lag")) %>%
  filter(!is.na(Rt))

reg_res4 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data4)

fitted_vals4 <- exp(fitted(reg_res4))
fitted_df4 <- reg_data4 %>%
  filter(!is.na(Rt)) %>%
  mutate(Rt_fitted = fitted_vals4)

dataset1_Simulx_new4_fitted <- dataset1_Simulx_new4 %>%
  left_join(fitted_df4 %>% rename(Rt_EE = Rt), by = c("dept_id", "day"))%>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94), labels = paste("Region", 1:94)))

# summary SEIRAHD 4 dataset with first days cut off
reg_data_m10 <- Rt_comp_res4$Rt_comp %>%
  mutate(lockdown1 = ifelse(between(day, 45, 99), 1, 0), 
         BG1 = ifelse(day > 99, 1, 0)) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         across(lockdown1:BG1, 
                function(x) lag(x, 5, default = 0), 
                .names = "{.col}_lag"))  %>%
  filter(day > 13)

reg_res_m10 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data_m10)

fitted_vals_m10 <- exp(fitted(reg_res_m10))
fitted_df_m10 <- reg_data_m10 %>%
  filter(day > 13) %>%
  mutate(Rt_fitted = fitted_vals_m10)

dataset1_Simulx_new4_m10_fitted <- dataset1_Simulx_new4 %>%
  left_join(fitted_df_m10 %>% rename(Rt_EE = Rt), by = c("dept_id", "day")) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94), labels = paste("Region", 1:94)))



# summary SEIRAHD 45dataset
reg_data5 <- Rt_comp_res5$Rt_comp %>%
  mutate(lockdown1 = ifelse(between(day, 45, 99), 1, 0), 
         BG1 = ifelse(day > 99, 1, 0))  %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         across(lockdown1:BG1, 
                function(x) lag(x, 5, default = 0), 
                .names = "{.col}_lag")) %>%
  filter(!is.na(Rt))

reg_res5 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data5)

fitted_vals5 <- exp(fitted(reg_res5))
fitted_df5 <- reg_data5 %>%
  filter(!is.na(Rt)) %>%
  mutate(Rt_fitted = fitted_vals5)

dataset1_Simulx_new5_fitted <- dataset1_Simulx_new5 %>%
  left_join(fitted_df5 %>% rename(Rt_EE = Rt), by = c("dept_id", "day"))%>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94), labels = paste("Region", 1:94)))

# summary SEIRAHD 4 dataset with first days cut off
reg_data5_m10 <- Rt_comp_res5$Rt_comp %>%
  mutate(lockdown1 = ifelse(between(day, 45, 99), 1, 0), 
         BG1 = ifelse(day > 99, 1, 0)) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         across(lockdown1:BG1, 
                function(x) lag(x, 5, default = 0), 
                .names = "{.col}_lag"))  %>%
  filter(day > 13)

reg_res5_m10 <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data5_m10)

fitted_vals5_m10 <- exp(fitted(reg_res5_m10))
fitted_df5_m10 <- reg_data5_m10 %>%
  filter(day > 13) %>%
  mutate(Rt_fitted = fitted_vals5_m10)

dataset1_Simulx_new5_m10_fitted <- dataset1_Simulx_new5 %>%
  left_join(fitted_df5_m10 %>% rename(Rt_EE = Rt), by = c("dept_id", "day")) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)), 
         dept_id2 = factor(dept_id2, levels = paste("Region", 1:94), labels = paste("Region", 1:94)))




## plots fits ---------------------------------------------------

rect_cols <- sequential_hcl(5, palette = "BluYl")

ggplot(dataset1_Simulx_new2_fitted %>% filter(dept_id %in% selected_depts & !is.na(dept_id2)), 
       aes(x = day, y = Rt_fitted)) + 
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = 43, y = Inf, label = "NPI 1", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  annotate("label", x = 96, y = Inf, label =  "NPI 2", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  geom_line(aes(linetype = "Regression fit Rt")) +
  geom_line(aes(y = Rt, linetype = "Real Rt")) +
  geom_line(aes(y = Rt_EE, linetype = "EpiEstim Rt")) +
  geom_xsideline(aes(y = IncI_unscaled)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4))+ 
  labs(linetype = "") + 
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 121)) + 
  labs(y = expression(R[t]), title = "Simulx infections, short NPI-free period, no days cut off, no NPI lag") +
  facet_wrap(~dept_id2) + 
  theme_bw() + 
  theme(legend.text = element_text(family = "serif", size = 13, hjust = 0), 
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))



ggplot(dataset1_Simulx_new3_fitted %>% filter(dept_id %in% selected_depts[1:9] & !is.na(dept_id2)), 
       aes(x = day)) + 
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = 43, y = Inf, label = "NPI 1", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  annotate("label", x = 96, y = Inf, label =  "NPI 2", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  geom_line(aes(y = Rt_fitted, linetype = "Regression fit Rt")) +
  geom_line(aes(y = Rt, linetype = "Real Rt")) +
  geom_line(aes(y = Rt_EE, linetype = "EpiEstim Rt")) +
  geom_xsideline(aes(y = IncI_unscaled)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4))+ 
  labs(linetype = "") + 
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 121)) + 
  labs(#title = "Simulx infections 3, short NPI-free period, no days cut off, no NPI lag", 
       y = expression(R[t]), x = "Day") +
  facet_wrap(~dept_id2) + 
  theme_bw() + 
  theme(legend.text = element_text(family = "serif", size = 13, hjust = 0), 
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))


ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/reg_fits_Simulx3_9reg_wo_title.jpeg", 
       dpi = 400, width = 10, height = 7.5)


ggplot(dataset1_Simulx_new4_fitted %>% filter(dept_id %in% selected_depts[1:9] & !is.na(dept_id2)), 
       aes(x = day, y = Rt_fitted)) + 
  annotate("rect", xmin = 45, xmax = 100, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 100, xmax = 151, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = 75, y = Inf, label = "NPI 1", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  annotate("label", x = 125, y = Inf, label =  "NPI 2", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  geom_line(aes(linetype = "Regression fit Rt")) +
  geom_line(aes(y = Rt, linetype = "Real Rt")) +
  geom_line(aes(y = Rt_EE, linetype = "EpiEstim Rt")) +
  geom_xsideline(aes(y = IncI_unscaled)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4))+ 
  labs(linetype = "") + 
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 151)) + 
  labs(#title = "Simulx infections, long NPI-free period, no days cut off, no NPI lag", 
       y = expression(R[t]), x = "Day") +
  facet_wrap(~dept_id2) + 
  theme_bw() + 
  theme(legend.text = element_text(family = "serif", size = 13, hjust = 0), 
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))

ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/reg_fits_Simulx4_9reg_wo_title.jpeg", 
       dpi = 400, width = 10, height = 7.5)



ggplot(dataset1_Simulx_new4_m10_fitted %>% filter(dept_id %in% selected_depts[1:9] & !is.na(dept_id2)),
       aes(x = day, y = Rt_fitted)) + 
  annotate("rect", xmin = 45, xmax = 100, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 100, xmax = 151, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = 75, y = Inf, label = "NPI 1", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  annotate("label", x = 125, y = Inf, label =  "NPI 2", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  geom_line(aes(linetype = "Regression fit Rt")) +
  geom_line(aes(y = Rt, linetype = "Real Rt")) +
  geom_line(aes(y = Rt_EE, linetype = "EpiEstim Rt")) +
  geom_xsideline(aes(y = IncI_unscaled)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4))+ 
  labs(linetype = "") + 
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 151)) + 
  labs(#title = "Simulx infections, long NPI-free period, first 10 days cut off, no NPI lag", 
       y = expression(R[t]), x = "Day") +
  facet_wrap(~dept_id2) + 
  theme_bw() + 
  theme(legend.text = element_text(family = "serif", size = 13, hjust = 0), 
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))

ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/reg_fits_Simulx4_m10_9reg_wo_title.jpeg", 
       dpi = 400, width = 10, height = 7.5)


ggplot(dataset1_Simulx_new5_fitted %>% filter(dept_id %in% selected_depts[1:9] & !is.na(dept_id2)), 
       aes(x = day, y = Rt_fitted)) + 
  annotate("rect", xmin = 45, xmax = 100, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 100, xmax = 151, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = 75, y = Inf, label = "NPI 1", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  annotate("label", x = 125, y = Inf, label =  "NPI 2", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  geom_line(aes(linetype = "Regression fit Rt")) +
  geom_line(aes(y = Rt, linetype = "Real Rt")) +
  geom_line(aes(y = Rt_EE, linetype = "EpiEstim Rt")) +
  geom_xsideline(aes(y = IncI_unscaled)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4))+ 
  labs(linetype = "") + 
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 151)) + 
  labs(#title = "Simulx infections, long NPI-free period, no days cut off, no NPI lag", 
    y = expression(R[t]), x = "Day") +
  facet_wrap(~dept_id2) + 
  theme_bw() + 
  theme(legend.text = element_text(family = "serif", size = 13, hjust = 0), 
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))

ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/reg_fits_Simulx5_9reg_wo_title.jpeg", 
       dpi = 400, width = 10, height = 7.5)



ggplot(dataset1_Simulx_new5_m10_fitted %>% filter(dept_id %in% selected_depts[1:9] & !is.na(dept_id2)),
       aes(x = day, y = Rt_fitted)) + 
  annotate("rect", xmin = 45, xmax = 100, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 100, xmax = 151, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = 75, y = Inf, label = "NPI 1", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  annotate("label", x = 125, y = Inf, label =  "NPI 2", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  geom_line(aes(linetype = "Regression fit Rt")) +
  geom_line(aes(y = Rt, linetype = "Real Rt")) +
  geom_line(aes(y = Rt_EE, linetype = "EpiEstim Rt")) +
  geom_xsideline(aes(y = IncI_unscaled)) + 
  ggside(scales = "free_y")  + 
  scale_xsidey_continuous(minor_breaks = NULL, breaks = scales::extended_breaks(n = 4))+ 
  labs(linetype = "") + 
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 151)) + 
  labs(#title = "Simulx infections, long NPI-free period, first 10 days cut off, no NPI lag", 
    y = expression(R[t]), x = "Day") +
  facet_wrap(~dept_id2) + 
  theme_bw() + 
  theme(legend.text = element_text(family = "serif", size = 13, hjust = 0), 
        ggside.panel.scale.x = .3, 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 13), 
        ggside.axis.text.y = element_text(family = "serif", size = 10))

ggsave("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots/reg_fits_Simulx5_m10_9reg_wo_title.jpeg", 
       dpi = 400, width = 10, height = 7.5)



## bias assessment  --------------------------------------------------------
bias_first_days2 <- fitted_df2 %>%
  group_by(dept_id) %>%
  filter(between(day, 8, 15)) %>%
  summarize(bias = mean(Rt_fitted - Rt_real))

bias_first_days3 <- fitted_df3 %>%
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

bias_first_days5 <- fitted_df5 %>%
  group_by(dept_id) %>%
  filter(between(day, 8, 44)) %>%
  summarize(bias = mean(Rt_fitted - Rt_real))

bias_first_days5_m10 <- fitted_df5_m10 %>%
  group_by(dept_id) %>%
  filter(between(day, 8, 44)) %>%
  summarize(bias = mean(Rt_fitted - Rt_real))

bias_first_days <- bias_first_days2 %>% mutate(model = "short") %>%
  bind_rows(bias_first_days3 %>% mutate(model = "short 3")) %>%
  bind_rows(bias_first_days4 %>% mutate(model = "long")) %>%
  bind_rows(bias_first_days_m10 %>% mutate(model = "long m10")) %>%
  bind_rows(bias_first_days5 %>% mutate(model = "long 3")) %>%
  bind_rows(bias_first_days5_m10 %>% mutate(model = "long 3 m10")) %>%
  mutate(model = as.factor(model))

ggplot(bias_first_days, aes(x = as.factor(dept_id), y = bias, fill = model)) + 
  geom_col(position = position_dodge(width = 0.8)) +
  scale_fill_brewer(palette = "Dark2")



## regression coefficients -------------------------------------------------
list_reg_res <- list(reg_res2, reg_res3, reg_res4, reg_res_m10, reg_res5, reg_res5_m10)
model_names <- c("Simulx 2", "Simulx 3", "Simulx 4", "Simulx 4 m10", 
                 "Simulx 5", "Simulx 5 m10")

reg_res_all <- list()
for(i in 1:6){
  coefs_Rt_reg <- coefficients(list_reg_res[[i]])$dept_id
  
  confint_Rt_reg <- data.frame(confint(list_reg_res[[i]], method="Wald"))[-c(1:2), ]
  names(confint_Rt_reg) <- c("CI_LL", "CI_UL")
  
  
  coefs_Rt_reg <- coefs_Rt_reg %>%
    dplyr::select(-1) %>%
    unique() %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
    cbind(confint_Rt_reg[-1,]) %>%
    mutate(parameter = factor(parameter,
                              levels = c("lockdown1", "BG1"),
                              labels = c("NPI 1", "NPI 2"))) 
  
  reg_res_all[[i]] <- coefs_Rt_reg %>%
    mutate(model = model_names[i])
}

reg_res_df <- do.call("rbind.data.frame", reg_res_all)
reg_res_df

ggplot(reg_res_df, aes(x = 1, y = value, ymin = CI_LL, ymax = CI_UL, col = model)) + 
  geom_pointrange(position = position_dodge(width = 1))  +
  geom_hline(data = data.frame(yval = -1.45, parameter = "NPI 1"), aes(yintercept = yval),
             linetype = "dashed", col = "black", linewidth = 0.8) +
  geom_hline(data = data.frame(yval = -0.5, parameter = "NPI 2"), aes(yintercept = yval),
             linetype = "dashed", col = "black", linewidth = 0.8) + 
  facet_wrap(~parameter) +
  labs(x = "", y = "coefficient value", 
       title = "Regression estimates") + 
  scale_color_brewer(palette = "Dark2") +
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.text = element_text(family = "serif", size = 13, hjust = 0),
        legend.title = element_text(family = "serif", size = 14), 
        plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13),  
        strip.text = element_text(family = "serif", size = 13), 
        axis.text.y = element_text(family = "serif", size = 12)) 


## bias metrics -------------------------------------------------
models_sel <- c("Simulx 3", "Simulx 5", "Simulx 5 m10")

metric_df_0.8 <- reg_res_df %>%
  filter(model %in% models_sel) %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.8), 
         CI_covers = ifelse(between(true_value, CI_LL, CI_UL), 1, 0), 
         bias = true_value - value, 
         rel_bias = abs(true_value - value)/abs(true_value)*100) %>%
  select(parameter, model, CI_covers, bias, rel_bias) %>%
  rename(coverage = CI_covers, relbias = rel_bias) %>%
  pivot_longer(cols = c(coverage, bias, relbias), 
               names_to = "metric", 
               values_to = "value") %>%
  pivot_wider(names_from = model, values_from = value) %>%
  mutate(metric = factor(metric, levels = c("coverage", "bias", "relbias"))) %>%
  arrange (parameter, metric) %>%
  mutate(metric = case_when(metric == "bias" ~ "absolute bias", 
                            metric == "coverage" ~ "CI coverage (%)", 
                            metric == "relbias" ~ "relative bias (%)"))

# print metrics table
metric_df_0.8 %>% 
  ungroup() %>%
  select(-parameter) %>%
  kable(digits = 2, format = "latex") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
  pack_rows("NPI 1", 1, 3) %>%
  pack_rows("NPI 2", 4, 6)

