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


# load data and calculate Rt ----------------------------------------------
popsize_df <- read.csv(paste0(dir2, "/popsize_df.csv")) %>%
  mutate(dept_id = ifelse(dept_id > 20, dept_id-1, dept_id))

indparams_Simulx2_1 <- read.table(paste0(dir2, "/ind_params/ind_2params_1.txt"), 
                                  header = TRUE, sep = " ")


data_list_Simulx_new2 <- loadRData(paste(dir2, "sim_res_Simulx_2params_new2_list.RData", sep = "/"))

dataset1_Simulx_new2 <- data_list_Simulx_new2[[1]] %>%
  left_join(popsize_df, by = c("id" = "dept_id")) %>%
  left_join(indparams_Simulx2_1, by = "id") %>%
  mutate(Rt = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
         IncE = transmission*S*(I+0.55*A)/10^4,
         IncE_unscaled = IncE*popsize/10^4,
         IncI_unscaled = IncI*popsize/10^4, 
         IncH_unscaled = IncH*popsize/10^4,
         lockdown1 = ifelse(between(time, 16, 70), 1, 0), 
         BG1 = ifelse(time > 70, 1, 0)) %>% 
  rename(dept_id = id, day = time)

true_Rt_df_Simulx_new2 <- dataset1_Simulx_new2 %>% dplyr::select(dept_id, day, Rt) %>% 
  rename(Rt_real = Rt)



# regressions with IncI and IncE
cl <- makeCluster(10)
registerDoParallel(cl)

Rt_comp_res_I <- Rt_calc_fun(Inc_df = dataset1_Simulx_new2, id_col_name = "dept_id", time_col_name = "day", 
                            Inc_col_name = "IncI_unscaled", model_name = "lag 0", 
                            Rt_ref_df = true_Rt_df_Simulx_new2, 
                            meansi = 10.1, stdsi = 8.75, Rt_prior = 1)

Rt_comp_res_E <- Rt_calc_fun(Inc_df = dataset1_Simulx_new2, id_col_name = "dept_id", time_col_name = "day", 
                            Inc_col_name = "IncE_unscaled", model_name = "lag 0", 
                            Rt_ref_df = true_Rt_df_Simulx_new2, 
                            meansi = 10.1, stdsi = 8.75, Rt_prior = 1)
stopCluster(cl)

reg_data_I <- Rt_comp_res_I$Rt_comp %>%
  mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
         BG1 = ifelse(day > 70, 1, 0)) %>%
  filter(!is.na(Rt))

reg_res_I <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data_I)

fitted_vals_I <- exp(fitted(reg_res_I))
fitted_df_I <- reg_data_I %>%
  filter(!is.na(Rt)) %>%
  mutate(Rt_fitted = fitted_vals_I)



reg_data_E <- Rt_comp_res_E$Rt_comp %>%
  mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
         BG1 = ifelse(day > 70, 1, 0)) %>%
  filter(!is.na(Rt))

reg_res_E <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data_E)

fitted_vals_E <- exp(fitted(reg_res_E))
fitted_df_E <- reg_data_E %>%
  filter(!is.na(Rt)) %>%
  mutate(Rt_fitted = fitted_vals_E)


coefs_I <- coefficients(reg_res_I)$dept_id[1, 2:3]
coefs_E <- coefficients(reg_res_E)$dept_id[1, 2:3]


confint_I <- data.frame(confint(reg_res_I, method="Wald"))[-c(1:2), ]
names(confint_I) <- c("CI_LL", "CI_UL")

confint_E <- data.frame(confint(reg_res_E, method="Wald"))[-c(1:2), ]
names(confint_E) <- c("CI_LL", "CI_UL")


reg_Rt_I <- coefs_I %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  cbind(confint_I[-1,]) %>%
  mutate(parameter = factor(parameter,
                            levels = c("lockdown1", "BG1"),
                            labels = c("NPI 1", "NPI 2"))) %>%
  mutate(model = "IncI")

reg_Rt_E <- coefs_E %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  cbind(confint_E[-1,]) %>%
  mutate(parameter = factor(parameter,
                            levels = c("lockdown1", "BG1"),
                            labels = c("NPI 1", "NPI 2"))) %>%
  mutate(model = "IncE")

reg_Rt_I
reg_Rt_E


# fitted values
fitted_vals_comp <- fitted_df_I %>% mutate(model = "IncI") %>%
  bind_rows(fitted_df_E %>% mutate(model = "IncE")) %>%
  mutate(dept_id2 = paste("Region", as.character(dept_id)))


set.seed(123)
selected_depts <- sample(1:94, size = 16, replace = FALSE)
rect_cols <- sequential_hcl(5, palette = "BluYl")

ggplot(fitted_vals_comp %>% filter(dept_id %in% selected_depts[1:9]), 
       aes(x = day, y = Rt_fitted, col = model)) + 
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = 43, y = Inf, label = "NPI 1", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  annotate("label", x = 96, y = Inf, label =  "NPI 2", 
           hjust = 0.5, vjust = 1, size = 4, fontface = 2, family = "serif") + 
  geom_line(aes(linetype = "Regression fit Rt"), linewidth = 0.8) +
  geom_line(aes(y = Rt_real, linetype = "Real Rt"), col = "black", linewidth = 0.8) +
  geom_line(aes(y = Rt, linetype = "EpiEstim Rt"), linewidth = 0.8) +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 121)) + 
  labs(y = expression(R[t]), title = "E vs. I as EpiEstim input", 
       linetype = "", col = "Rt est. input") +
  facet_wrap(~dept_id2) + 
  theme_bw()


# bias in regression values
metric_table_IncE <-  reg_Rt_I %>%
  bind_rows(reg_Rt_E) %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.5), 
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

metric_table_IncE %>% 
  select(-parameter) %>%
  kable(digits = 2, format = "html") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
  pack_rows("NPI 1", 1, 3) %>%
  pack_rows("NPI 2", 4, 6)

