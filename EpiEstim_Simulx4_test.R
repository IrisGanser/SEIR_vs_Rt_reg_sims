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
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/deSolve_SIR_model_function.R")
source("~/PhD/COVID_France/Dropbox_iris_covid/departement/Données_SPF/Data/data_functions.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")


dir2 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"


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
         IncI_unscaled = round(IncI*popsize/10^4), 
         IncH_unscaled = round(IncH*popsize/10^4), 
         lockdown1 = ifelse(between(time, 45, 99), 1, 0), 
         BG1 = ifelse(time > 99, 1, 0)) %>% 
  rename(dept_id = id, day = time)

true_Rt_df_Simulx_new4 <- dataset1_Simulx_new4 %>% dplyr::select(dept_id, day, Rt_SEIRAHD) %>% 
  rename(Rt_real = Rt_SEIRAHD)



ggplot(dataset1_Simulx_new4, aes(x = day, y = Rt_SEIRAHD, group = dept_id)) + 
  geom_line(aes(col = "SEIRAHD"), alpha = 0.5) +
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  labs(title = "Rt from Simulx 4 model", y = "Rt", col = "") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() 


# Rt function
Rt_calc_fun <- function(Inc_df, id_col_name, time_col_name, Inc_col_name, model_name, 
                        tstart = NULL, tend = NULL, Rt_ref_df, 
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
        mutate(t = ceiling((t_end+t_start)/2)-1,
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
        mutate(t = ceiling((t_end+t_start)/2),
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


Rt_reg_only_fun <- function(data_for_est, lag_NPIs = FALSE, lag_days = 0, model_name, 
                            fits = FALSE){
  
  if(lag_NPIs){
    
    reg_data <- data_for_est %>%
      dplyr::select(dept_id, day, Rt, lockdown1, BG1) %>%
      group_by(dept_id) %>%
      mutate(lockdown1 = lag(lockdown1, n = lag_days, default = 0), 
             BG1 = lag(BG1, n = lag_days, default = 0)) %>%
      ungroup()
    
    reg_res <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data)
    
  } else{
    
    reg_data <- data_for_est %>%
      dplyr::select(dept_id, day, Rt, lockdown1, BG1) 
    
    reg_res <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data)
    
  }
  
  
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
                              labels = c("NPI 1", "NPI 2"))) %>%
    mutate(lag = lag_days, model = model_name)
  
  
  if(fits){
    fitted_vals <- exp(fitted(reg_res))
    fitted_df <- data_for_est %>%
      filter(!is.na(Rt)) %>%
      mutate(Rt_fitted = fitted_vals, lag = lag_days, model = model_name)
    
    return(list(coefs = coefs_Rt_reg, fits = fitted_df))
  } else{
    
    return(coefs_Rt_reg)
  }
  
}



#### NPI lag ####
lag <- 0:15
new_tstart <- list(2:150, 2:149, 2:148, 2:147, 2:146, 2:145, 2:144)
new_tend <- list(2:150, 3:150, 4:150, 5:150, 6:150, 7:150, 8:150)
interval <- 1:7

cl <- makeCluster(6)
registerDoParallel(cl)

# Simulx
list_Rt_reg_lag_Simulx <- list()
list_Rt_reg_lag_fits_Simulx <- list()
list_Rt_reg_lag_Simulx_m10 <- list()
list_Rt_reg_lag_fits_Simulx_m10 <- list()

for(i in 1:length(lag)){
  Rt_comp_res <- Rt_calc_fun(Inc_df = dataset1_Simulx_new4, id_col_name = "dept_id", time_col_name = "day", 
                             Inc_col_name = "IncI_unscaled", model_name = lag[i], 
                             Rt_ref_df = true_Rt_df_Simulx_new4, tstart = new_tstart[[7]], tend = new_tend[[7]], 
                             meansi = 10.1, stdsi = 8.75, Rt_prior = 1)
  
  reg <- Rt_reg_only_fun(data_for_est = Rt_comp_res$Rt_comp %>%
                           mutate(lockdown1 = ifelse(between(day, 45, 99), 1, 0), 
                                  BG1 = ifelse(day > 99, 1, 0)), 
                         model_name = lag[i], lag_NPIs = TRUE, lag_days = lag[i], fits = TRUE)
  
  reg_m10 <- Rt_reg_only_fun(data_for_est = Rt_comp_res$Rt_comp %>%
                               filter(day > 10) %>%
                               mutate(lockdown1 = ifelse(between(day, 45, 99), 1, 0), 
                                      BG1 = ifelse(day > 99, 1, 0)), 
                             model_name = lag[i], lag_NPIs = TRUE, lag_days = lag[i], fits = TRUE)
  
  list_Rt_reg_lag_Simulx[[i]] <- reg$coefs
  list_Rt_reg_lag_fits_Simulx[[i]] <- reg$fits
  
  list_Rt_reg_lag_Simulx_m10[[i]] <- reg_m10$coefs
  list_Rt_reg_lag_fits_Simulx_m10[[i]] <- reg_m10$fits
  
}

reg_res_lag_Simulx_df <- do.call("rbind.data.frame", list_Rt_reg_lag_Simulx)
reg_res_lag_Simulx_fits_df <- do.call("rbind.data.frame", list_Rt_reg_lag_fits_Simulx)
reg_res_lag_Simulx_m10_df <- do.call("rbind.data.frame", list_Rt_reg_lag_Simulx_m10)
reg_res_lag_Simulx_fits_m10_df <- do.call("rbind.data.frame", list_Rt_reg_lag_fits_Simulx_m10)


list_Rt_reg_lag_Simulx_h <- list()
list_Rt_reg_lag_fits_Simulx_h <- list()
list_Rt_reg_lag_Simulx_h_m10 <- list()
list_Rt_reg_lag_fits_Simulx_h_m10 <- list()

for(i in 1:length(lag)){
  Rt_comp_res <- Rt_calc_fun(Inc_df = dataset1_Simulx_new4, id_col_name = "dept_id", time_col_name = "day", 
                             Inc_col_name = "IncH_unscaled", model_name = lag[i], 
                             Rt_ref_df = true_Rt_df_Simulx_new4, tstart = new_tstart[[7]], tend = new_tend[[7]], 
                             meansi = 10.1, stdsi = 8.75, Rt_prior = 1)
  
  reg <- Rt_reg_only_fun(data_for_est = Rt_comp_res$Rt_comp %>%
                           mutate(lockdown1 = ifelse(between(day, 45, 99), 1, 0), 
                                  BG1 = ifelse(day > 99, 1, 0)), 
                         model_name = lag[i], lag_NPIs = TRUE, lag_days = lag[i], fits = TRUE)
  
  
  reg_m10 <- Rt_reg_only_fun(data_for_est = Rt_comp_res$Rt_comp %>%
                               filter(day > 10) %>%
                           mutate(lockdown1 = ifelse(between(day, 45, 99), 1, 0), 
                                  BG1 = ifelse(day > 99, 1, 0)), 
                         model_name = lag[i], lag_NPIs = TRUE, lag_days = lag[i], fits = TRUE)
  
  
  list_Rt_reg_lag_Simulx_h[[i]] <- reg$coefs
  list_Rt_reg_lag_fits_Simulx_h[[i]] <- reg$fits
  list_Rt_reg_lag_Simulx_h_m10[[i]] <- reg_m10$coefs
  list_Rt_reg_lag_fits_Simulx_h_m10[[i]] <- reg_m10$fits
  
}

reg_res_lag_Simulx_h_df <- do.call("rbind.data.frame", list_Rt_reg_lag_Simulx_h)
reg_res_lag_Simulx_h_fits_df <- do.call("rbind.data.frame", list_Rt_reg_lag_fits_Simulx_h)
reg_res_lag_Simulx_h_m10_df <- do.call("rbind.data.frame", list_Rt_reg_lag_Simulx_h_m10)
reg_res_lag_Simulx_h_fits_m10_df <- do.call("rbind.data.frame", list_Rt_reg_lag_fits_Simulx_h_m10)

stopCluster(cl)


comp_df_lag_fits <- reg_res_lag_Simulx_fits_df %>%
  mutate(model2 = "Simulx IncI") %>%
  bind_rows(reg_res_lag_Simulx_fits_m10_df %>% mutate(model2 = "Simulx IncI minus 10")) %>%
  bind_rows(reg_res_lag_Simulx_h_fits_df %>% mutate(model2 = "Simulx IncH"))  %>%
  bind_rows(reg_res_lag_Simulx_h_fits_m10_df %>% mutate(model2 = "Simulx IncH minus 10")) 

reg_comp_df_lag <- reg_res_lag_Simulx_df %>%
  mutate(model2 = "Simulx IncI") %>%
  bind_rows(reg_res_lag_Simulx_m10_df %>% mutate(model2 = "Simulx IncI minus 10")) %>%
  bind_rows(reg_res_lag_Simulx_h_df %>% mutate(model2 = "Simulx IncH")) %>%
  bind_rows(reg_res_lag_Simulx_h_m10_df %>% mutate(model2 = "Simulx IncH minus 10")) %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.5))



ggplot(reg_comp_df_lag, aes(x = model, y = value, col = model2)) + 
  geom_pointrange(aes(ymin = CI_LL, ymax = CI_UL), position = position_dodge(width = 0.3)) + 
  geom_line(aes(y = true_value), linetype = "dashed", col = "darkred") +
  scale_x_continuous(breaks = 0:15) + 
  facet_wrap(~parameter) + 
  labs(y = "coefficient value", x = "NPI lag [days]", col = "model") + 
  scale_color_brewer(palette = "Paired") +
  theme_bw()


ggplot(comp_df_lag_fits %>% filter(dept_id %in% c(1, 4, 9, 13)), 
       aes(x = day, col = as.factor(lag))) + 
  geom_line(aes(y = Rt, linetype = "EpiEstim Rt"), col = "black", linewidth = 0.8) +
  geom_line(aes(y = Rt_real, linetype = "Real Rt"), col = "black", linewidth = 0.8) +
  geom_line(aes(y = Rt_fitted, linetype = "Regression fit Rt")) +
  facet_grid(cols = vars(dept_id), rows = vars(model2)) + 
  labs(col = "NPI lag [days]", linetype = "Rt") + 
  scale_color_viridis_d() +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) +
  theme_bw() 

ggsave("EpiEstim_test_Simulx4_reg_fits.jpeg", width = 14, height = 8)

