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

setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/Rt_trajectories")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/deSolve_SIR_model_function.R")
source("~/PhD/COVID_France/Dropbox_iris_covid/departement/Données_SPF/Data/data_functions.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")


dir2 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"
dir3 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once7"

popsize_df <- read.csv(paste0(dir2, "/popsize_df.csv")) %>%
  mutate(dept_id = ifelse(dept_id > 20, dept_id-1, dept_id))

indparams_Simulx4_1 <- read.table(paste0(dir2, "/ind_params/ind_2params_new4_1.txt"), 
                                  header = TRUE, sep = " ")


#### load data and calculate Rt ####
data_list_Simulx_new4 <- loadRData(paste(dir2, "sim_res_Simulx_2params_new4_list.RData", sep = "/"))

dataset1_Simulx_new4 <- data_list_Simulx_new4[[1]] %>%
  left_join(popsize_df, by = c("id" = "dept_id")) %>%
  left_join(indparams_Simulx4_1, by = "id") %>%
  mutate(Rt_SEIRAHD = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
         IncI_unscaled = IncI*popsize/10^4, 
         IncH_unscaled = IncH*popsize/10^4,
         lockdown1 = ifelse(between(time, 45, 99), 1, 0), 
         BG1 = ifelse(time > 99, 1, 0)) %>% 
  rename(dept_id = id, day = time)

true_Rt_df_Simulx_new4 <- dataset1_Simulx_new4 %>% dplyr::select(dept_id, day, Rt_SEIRAHD) %>% 
  rename(Rt_real = Rt_SEIRAHD)



data_ABM_rm_cov_long <- read.csv(paste0(dir3, "/data_covasim_rm7_Rt_1.csv")) 
data_ABM_hybrid_cov_long <- read.csv(paste0(dir3, "/data_covasim_hybrid7_Rt_1.csv"))

true_Rt_df_ABM_rm_long <- data_ABM_rm_cov_long %>% 
  dplyr::select(dept_id, day, Rt) %>% 
  rename(Rt_real = Rt)
true_Rt_df_ABM_hybrid_long <- data_ABM_hybrid_cov_long %>% 
  dplyr::select(dept_id, day, Rt) %>% 
  rename(Rt_real = Rt)


ggplot(dataset1_Simulx_new4, aes(x = day, y = Rt_SEIRAHD, group = dept_id)) + 
  geom_line(aes(col = "SEIRAHD"), alpha = 0.5) +
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  labs(title = "Rt from Simulx 4 model", y = "Rt", col = "") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() 


# Rt function
Rt_calc_fun <- function(Inc_df, id_col_name, time_col_name, Inc_col_name, model_name, 
                        tstart = NULL, tend = NULL, Rt_ref_df, cut_days = 0,
                        Rt_prior = 1, Rt_sd_prior = 2, meansi = 10.1, stdsi = 8.75){
  # select input data according to desired input for Rt calculation
  Inc_data <- Inc_df[c(id_col_name, time_col_name, Inc_col_name)]
  colnames(Inc_data) <- c("id", "time", "Inc")
  
  # Rt calculation
  Rt_list <- vector(mode = "list")
  
  if(is_null(tstart) | (is_null(tend))){
    Rt_list <- foreach(i = 1:94, .packages = c("EpiEstim", "tidyverse")) %dopar% {
      Rt_estim_I <- estimate_R(Inc_data$Inc[Inc_data$id == unique(Inc_data$id)[i]],
                               method = "parametric_si",
                               config = make_config(list(mean_si = meansi, 
                                                         std_si = stdsi, 
                                                         mean_prior = Rt_prior, 
                                                         std_prior = Rt_sd_prior)))$R
      Rt_estim_I <- Rt_estim_I %>%
        mutate(t = (t_end+t_start)/2 + cut_days,
               id = unique(Inc_data$id)[i], .before = t_start) %>%
        select(t, t_start, t_end, id, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`) %>%
        rename(Rt = `Mean(R)`, CI_LL = `Quantile.0.025(R)`, CI_UL = `Quantile.0.975(R)`)
    } 
  } 
  
  if(!is_null(tstart) & !(is_null(tend))){
    Rt_list <- foreach(i = 1:94, .packages = c("EpiEstim", "tidyverse")) %dopar% {
      Rt_estim_I <- estimate_R(Inc_data$Inc[Inc_data$id == unique(Inc_data$id)[i]],
                               method = "parametric_si",
                               config = make_config(list(
                                 mean_si = meansi, std_si = stdsi, 
                                 t_start = tstart, t_end = tend, 
                                 mean_prior = Rt_prior, 
                                 std_prior = Rt_sd_prior)))$R
      Rt_estim_I <- Rt_estim_I %>%
        mutate(t = (t_end+t_start)/2,
               id = unique(Inc_data$id)[i], .before = t_start) %>%
        select(t, t_start, t_end, id, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`) %>%
        rename(Rt = `Mean(R)`, CI_LL = `Quantile.0.025(R)`, CI_UL = `Quantile.0.975(R)`)
    } 
  }
  
  Rt_df <- do.call("rbind.data.frame", Rt_list) %>%
    rename(dept_id = id, day = t)
  
  Rt_comp_df <- Rt_df %>%
    left_join(., Rt_ref_df, by = c("dept_id", "day")) %>%
    mutate(Rt_residual = Rt - Rt_real, 
           model = model_name)
  
  RMSE_df <- Rt_comp_df %>%
    filter(!is.na(Rt)) %>%
    group_by(dept_id, model) %>%
    summarize(RMSE = sqrt(sum(Rt_residual^2)/n())) %>%
    ungroup() %>%
    mutate(RMSE_total = mean(RMSE))
  
  return(list(Rt_comp = Rt_comp_df, RMSE = RMSE_df))
  
}


#### regression tests ####
cl <- makeCluster(6)
registerDoParallel(cl)

Rt_comp_res <- Rt_calc_fun(Inc_df = dataset1_Simulx_new4, id_col_name = "dept_id", time_col_name = "day", 
                           Inc_col_name = "IncI_unscaled", model_name = "lag 0", 
                           Rt_ref_df = true_Rt_df_Simulx_new4, 
                           meansi = 10.1, stdsi = 8.75, Rt_prior = 1)
stopCluster(cl)

reg_data <- Rt_comp_res$Rt_comp %>%
  mutate(lockdown1 = ifelse(between(day, 45, 99), 1, 0), 
         BG1 = ifelse(day > 99, 1, 0)) %>%
  filter(!is.na(Rt))

reg_res <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data)

summary(reg_res)
coefs_reg_res <- coef(reg_res)$dept_id
hist(exp(coefs_reg_res[, 1]), breaks = 15)
range(coefs_reg_res[, 1])


fitted_vals <- exp(fitted(reg_res))
fitted_df <- reg_data %>%
  filter(!is.na(Rt)) %>%
  mutate(Rt_fitted = fitted_vals)

rect_cols <- sequential_hcl(5, palette = "BluYl")
set.seed(123)
selected_depts <- sample(1:94, size = 16, replace = FALSE)

ggplot(fitted_df %>% filter(dept_id %in% selected_depts), aes(x = day, y = Rt_fitted)) + 
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
  facet_wrap(~dept_id) + 
  theme_bw()

bias_first_days <- fitted_df %>%
  group_by(dept_id) %>%
  filter(between(day, 8, 44)) %>%
  summarize(bias = mean(Rt_fitted - Rt_real))

ggplot(bias_first_days, aes(x = dept_id, y = bias)) + 
  geom_col()

ggplot(bias_first_days, aes(x = bias)) + 
  geom_histogram(col = "black", bins = 20) + 
  labs(title = "Histogram of bias in random intercept estimation") + 
  theme_bw()

mean(bias_first_days$bias)
median(bias_first_days$bias)

ggplot(fitted_df %>% filter(day == 30), aes(x = Rt_real)) + 
  geom_histogram(aes(fill = "Real Rt"), alpha = 0.5, bins = 20, col = "black") + 
  geom_histogram(aes(x = Rt_fitted, fill = "Regression fit Rt"), alpha = 0.5, bins = 20, col = "black") +
  geom_histogram(aes(x = Rt, fill = "EpiEstim Rt"), alpha = 0.5, bins = 20, col = "black") +
  labs(x = "R0", title = "Comparison of real vs. fitted R0") +
  theme_bw()


fitted_df_long <- fitted_df %>%
  pivot_longer(cols = c(Rt_real, Rt, Rt_fitted), names_to = "Rt_type", values_to = "Rt") %>%
  mutate(Rt_type = case_when(Rt_type == "Rt" ~ "EpiEstim", 
                             Rt_type == "Rt_fitted" ~ "Regression fit", 
                             Rt_type == "Rt_real" ~ "Real"), 
         Rt_type = factor(Rt_type, levels = c("Real", "EpiEstim", "Regression fit")))


ggplot(fitted_df_long %>% filter(day == 30), aes(x = Rt, y = Rt_type, fill = Rt_type)) + 
  geom_density_ridges2(stat = "binline", bins = 20, scale = 0.95, draw_baseline = TRUE, alpha = 0.6) +
  labs(x = expression(R[0]), title = expression("Distribution of estimated vs. estimated"~ R[0]~"across departments"), 
       y = "", fill = expression(R[0]~"type")) +
  scale_y_discrete(limits = rev(levels(fitted_df_long$Rt_type))) + 
  scale_fill_brewer(palette = "Greens") + 
  theme_bw()
