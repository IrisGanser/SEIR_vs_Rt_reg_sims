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
dir3 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once7"

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

true_Rt_df_Simulx <- dataset1_Simulx %>% dplyr::select(dept_id, day, Rt_SEIRAHD) %>% 
  rename(Rt_real = Rt_SEIRAHD)


data_ABM_rm_cov <- read.csv(paste0(dir3, "/data_covasim_rm7_Rt_1.csv")) %>%
  filter(day > 29) %>%
  mutate(day = day - 29)
data_ABM_hybrid_cov <- read.csv(paste0(dir3, "/data_covasim_hybrid7_Rt_1.csv")) %>%
  filter(day > 29) %>%
  mutate(day = day - 29)

true_Rt_df_ABM_rm <- data_ABM_rm_cov %>% dplyr::select(dept_id, day, Rt) %>% 
  rename(Rt_real = Rt)
true_Rt_df_ABM_hybrid <- data_ABM_hybrid_cov %>% dplyr::select(dept_id, day, Rt) %>% 
  rename(Rt_real = Rt)


ggplot(dataset1_Simulx, aes(x = day, y = Rt_SEIRAHD, group = dept_id)) + 
  geom_line(aes(col = "SEIRAHD"), alpha = 0.5) +
  geom_line(data = data_ABM_rm_cov, aes(col = "ABM rm", y = Rt), alpha = 0.5) +
  geom_line(data = data_ABM_hybrid_cov, aes(col = "ABM hybrid", y = Rt), alpha = 0.5) +
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  labs(title = "Rt from three data generation models", y = "Rt", col = "") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() 


# Rt function
Rt_calc_fun <- function(Inc_df, id_col_name, time_col_name, Inc_col_name, model_name, 
                        tstart = NULL, tend = NULL, Rt_ref_df, 
                        Rt_prior = 2, Rt_sd_prior = 3, meansi = 10.1, stdsi = 8.75){
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


Rt_calc_ad_fun <- function(Inc_df, id_col_name, time_col_name, Inc_col_name, model_name, 
                        tstart = NULL, tend = NULL, Rt_ref_df, 
                        Rt_prior = 2, Rt_sd_prior = 3, meansi = 10.1, stdsi = 8.75, 
                        day_assigned){
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
        mutate(t = t_start + day_assigned,
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
        mutate(t = t_start + day_assigned,
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


  
#### different smoothing window lengths ####
cl <- makeCluster(8)
registerDoParallel(cl)

new_tstart <- list(2:121, 2:120, 2:119, 2:118, 2:117, 2:116, 2:115)
new_tend <- list(2:121, 3:121, 4:121, 5:121, 6:121, 7:121, 8:121)
interval <- 1:7

# Simulx
list_Rt_res_sm_Simulx <- list()
list_RMSE_sm_Simulx <- list()
list_Rt_reg_sm_Simulx <- list()
for(i in 1:length(interval)){
  Rt_comp_res <- Rt_calc_fun(Inc_df = dataset1_Simulx, id_col_name = "dept_id", time_col_name = "day", 
                             Inc_col_name = "IncI_unscaled", model_name = interval[i], 
                             Rt_ref_df = true_Rt_df_Simulx, tstart = new_tstart[[i]], tend = new_tend[[i]], 
                             meansi = 10.1, stdsi = 8.75, Rt_prior = 1, Rt_sd_prior = 2)
  
  
  list_RMSE_sm_Simulx[[i]] <- Rt_comp_res$RMSE
  list_Rt_res_sm_Simulx[[i]] <- Rt_comp_res$Rt_comp 
  
  
  list_Rt_reg_sm_Simulx[[i]] <- Rt_reg_only_fun(data_for_est = Rt_comp_res$Rt_comp %>%
                                             mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
                                                    BG1 = ifelse(day > 70, 1, 0)), 
                                           model_name = interval[i])
}
  
Rt_res_sm_Simulx_df <- do.call("rbind.data.frame", list_Rt_res_sm_Simulx)
reg_res_sm_Simulx_df <- do.call("rbind.data.frame", list_Rt_reg_sm_Simulx)


# ABM rm
list_Rt_res_sm_ABM_rm <- list()
list_RMSE_sm_ABM_rm <- list()
list_Rt_reg_sm_ABM_rm <- list()
for(i in 1:length(interval)){
  Rt_comp_res <- Rt_calc_fun(Inc_df = data_ABM_rm_cov, id_col_name = "dept_id", time_col_name = "day", 
                             Inc_col_name = "IncI", model_name = interval[i], 
                             Rt_ref_df = true_Rt_df_ABM_rm, tstart = new_tstart[[i]], tend = new_tend[[i]], 
                             meansi = 8.2, stdsi = 5, Rt_prior = 1, Rt_sd_prior = 2)
  
  
  list_RMSE_sm_ABM_rm[[i]] <- Rt_comp_res$RMSE
  list_Rt_res_sm_ABM_rm[[i]] <- Rt_comp_res$Rt_comp 
  
  
  list_Rt_reg_sm_ABM_rm[[i]] <- Rt_reg_only_fun(data_for_est = Rt_comp_res$Rt_comp %>%
                                                  mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
                                                         BG1 = ifelse(day > 70, 1, 0)), 
                                                model_name = interval[i])
}

Rt_res_sm_ABM_rm_df <- do.call("rbind.data.frame", list_Rt_res_sm_ABM_rm)
reg_res_sm_ABM_rm_df <- do.call("rbind.data.frame", list_Rt_reg_sm_ABM_rm)


# ABM hybrid
list_Rt_res_sm_ABM_hybrid <- list()
list_RMSE_sm_ABM_hybrid <- list()
list_Rt_reg_sm_ABM_hybrid <- list()
for(i in 1:length(interval)){
  Rt_comp_res <- Rt_calc_fun(Inc_df = data_ABM_hybrid_cov, id_col_name = "dept_id", time_col_name = "day", 
                             Inc_col_name = "IncI", model_name = interval[i], 
                             Rt_ref_df = true_Rt_df_ABM_hybrid, tstart = new_tstart[[i]], tend = new_tend[[i]], 
                             meansi = 8.2, stdsi = 5, Rt_prior = 1, Rt_sd_prior = 2)
  
  
  list_RMSE_sm_ABM_hybrid[[i]] <- Rt_comp_res$RMSE
  list_Rt_res_sm_ABM_hybrid[[i]] <- Rt_comp_res$Rt_comp 
  
  
  list_Rt_reg_sm_ABM_hybrid[[i]] <- Rt_reg_only_fun(data_for_est = Rt_comp_res$Rt_comp %>%
                                                  mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
                                                         BG1 = ifelse(day > 70, 1, 0)), 
                                                model_name = interval[i])
}

Rt_res_sm_ABM_hybrid_df <- do.call("rbind.data.frame", list_Rt_res_sm_ABM_hybrid)
reg_res_sm_ABM_hybrid_df <- do.call("rbind.data.frame", list_Rt_reg_sm_ABM_hybrid)

  
stopCluster(cl)
  

comp_df_sm <- Rt_res_sm_Simulx_df %>%
  mutate(model2 = "Simulx") %>%
  bind_rows(Rt_res_sm_ABM_hybrid_df %>% mutate(model2 = "ABM hybrid")) %>%
  bind_rows(Rt_res_sm_ABM_rm_df %>% mutate(model2 = "ABM rm")) 

reg_comp_df_sm <- reg_res_sm_Simulx_df %>%
  mutate(model2 = "Simulx") %>%
  bind_rows(reg_res_sm_ABM_hybrid_df %>% mutate(model2 = "ABM hybrid")) %>%
  bind_rows(reg_res_sm_ABM_rm_df %>% mutate(model2 = "ABM rm")) %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.5))


ggplot(reg_comp_df_sm, aes(x = model, y = value, col = model2)) + 
  geom_pointrange(aes(ymin = CI_LL, ymax = CI_UL), position = position_dodge(width = 0.3)) + 
  geom_line(aes(y = true_value), linetype = "dashed", col = "darkred") +
  scale_x_continuous(breaks = 1:7) + 
  facet_wrap(~parameter) + 
  labs(y = "coefficient value", x = "Smoothing window", col = "model") + 
  scale_color_brewer(palette = "Dark2") +
  theme_bw()

ggplot(comp_df_sm %>% filter(dept_id %in% c(1, 4, 9, 13)), 
       aes(x = day, y = Rt)) + 
  geom_ribbon(aes(ymin = CI_LL, ymax = CI_UL, fill = as.factor(model)), alpha = 0.1) +
  geom_line(aes(linetype = "estimated", col = as.factor(model))) + 
  geom_line(aes(y = Rt_real, linetype = "real")) + 
  facet_grid(cols = vars(dept_id), rows = vars(model2)) + 
  scale_color_viridis_d() + 
  scale_fill_viridis_d() +
  labs(linetype = "", col = "smoothing window", fill = "smoothing window") +
  theme_bw()
  



#### assignment day ####
cl <- makeCluster(8)
registerDoParallel(cl)

assignment_day <- 0:6

# Simulx
list_Rt_res_ad_Simulx <- list()
list_RMSE_ad_Simulx <- list()
list_Rt_reg_ad_Simulx <- list()
for(i in 1:length(interval)){
  Rt_comp_res <- Rt_calc_ad_fun(Inc_df = dataset1_Simulx, id_col_name = "dept_id", time_col_name = "day", 
                             Inc_col_name = "IncI_unscaled", model_name = assignment_day[i] + 1, 
                             Rt_ref_df = true_Rt_df_Simulx, tstart = new_tstart[[7]], tend = new_tend[[7]], 
                             meansi = 10.1, stdsi = 8.75, Rt_prior = 1, Rt_sd_prior = 2, 
                             day_assigned = assignment_day[i])
  
  
  list_RMSE_ad_Simulx[[i]] <- Rt_comp_res$RMSE
  list_Rt_res_ad_Simulx[[i]] <- Rt_comp_res$Rt_comp 
  
  
  list_Rt_reg_ad_Simulx[[i]] <- Rt_reg_only_fun(data_for_est = Rt_comp_res$Rt_comp %>%
                                                  mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
                                                         BG1 = ifelse(day > 70, 1, 0)), 
                                                model_name = interval[i])
}

Rt_res_ad_Simulx_df <- do.call("rbind.data.frame", list_Rt_res_ad_Simulx)
reg_res_ad_Simulx_df <- do.call("rbind.data.frame", list_Rt_reg_ad_Simulx)


# ABM rm
list_Rt_res_ad_ABM_rm <- list()
list_RMSE_ad_ABM_rm <- list()
list_Rt_reg_ad_ABM_rm <- list()
for(i in 1:length(interval)){
  Rt_comp_res <- Rt_calc_ad_fun(Inc_df = data_ABM_rm_cov, id_col_name = "dept_id", time_col_name = "day", 
                             Inc_col_name = "IncI", model_name = assignment_day[i] + 1, 
                             Rt_ref_df = true_Rt_df_ABM_rm, tstart = new_tstart[[7]], tend = new_tend[[7]], 
                             meansi = 8.2, stdsi = 5, Rt_prior = 1, Rt_sd_prior = 2, 
                             day_assigned = assignment_day[i])
  
  
  list_RMSE_ad_ABM_rm[[i]] <- Rt_comp_res$RMSE
  list_Rt_res_ad_ABM_rm[[i]] <- Rt_comp_res$Rt_comp 
  
  
  list_Rt_reg_ad_ABM_rm[[i]] <- Rt_reg_only_fun(data_for_est = Rt_comp_res$Rt_comp %>%
                                                  mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
                                                         BG1 = ifelse(day > 70, 1, 0)), 
                                                model_name = interval[i])
}

Rt_res_ad_ABM_rm_df <- do.call("rbind.data.frame", list_Rt_res_ad_ABM_rm)
reg_res_ad_ABM_rm_df <- do.call("rbind.data.frame", list_Rt_reg_ad_ABM_rm)


# ABM hybrid
list_Rt_res_ad_ABM_hybrid <- list()
list_RMSE_ad_ABM_hybrid <- list()
list_Rt_reg_ad_ABM_hybrid <- list()
for(i in 1:length(interval)){
  Rt_comp_res <- Rt_calc_ad_fun(Inc_df = data_ABM_hybrid_cov, id_col_name = "dept_id", time_col_name = "day", 
                             Inc_col_name = "IncI", model_name = assignment_day[i] + 1, 
                             Rt_ref_df = true_Rt_df_ABM_hybrid, tstart = new_tstart[[7]], tend = new_tend[[7]], 
                             meansi = 8.2, stdsi = 5, Rt_prior = 1, Rt_sd_prior = 2, 
                             day_assigned = assignment_day[i])
  
  
  list_RMSE_ad_ABM_hybrid[[i]] <- Rt_comp_res$RMSE
  list_Rt_res_ad_ABM_hybrid[[i]] <- Rt_comp_res$Rt_comp 
  
  
  list_Rt_reg_ad_ABM_hybrid[[i]] <- Rt_reg_only_fun(data_for_est = Rt_comp_res$Rt_comp %>%
                                                      mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
                                                             BG1 = ifelse(day > 70, 1, 0)), 
                                                    model_name = interval[i])
}

Rt_res_ad_ABM_hybrid_df <- do.call("rbind.data.frame", list_Rt_res_ad_ABM_hybrid)
reg_res_ad_ABM_hybrid_df <- do.call("rbind.data.frame", list_Rt_reg_ad_ABM_hybrid)


stopCluster(cl)


comp_df_ad <- Rt_res_ad_Simulx_df %>%
  mutate(model2 = "Simulx") %>%
  bind_rows(Rt_res_ad_ABM_hybrid_df %>% mutate(model2 = "ABM hybrid")) %>%
  bind_rows(Rt_res_ad_ABM_rm_df %>% mutate(model2 = "ABM rm")) 

reg_comp_df_ad <- reg_res_ad_Simulx_df %>%
  mutate(model2 = "Simulx") %>%
  bind_rows(reg_res_ad_ABM_hybrid_df %>% mutate(model2 = "ABM hybrid")) %>%
  bind_rows(reg_res_ad_ABM_rm_df %>% mutate(model2 = "ABM rm")) %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.5))


ggplot(reg_comp_df_ad, aes(x = model, y = value, col = model2)) + 
  geom_pointrange(aes(ymin = CI_LL, ymax = CI_UL), position = position_dodge(width = 0.3)) + 
  geom_line(aes(y = true_value), linetype = "dashed", col = "darkred") +
  scale_x_continuous(breaks = 1:7) + 
  facet_wrap(~parameter) + 
  labs(y = "coefficient value", x = "Rt assignment day", col = "model") + 
  scale_color_brewer(palette = "Dark2") +
  theme_bw()

ggplot(comp_df_ad %>% filter(dept_id %in% c(1, 4, 9, 13)), 
       aes(x = day, y = Rt)) + 
  geom_ribbon(aes(ymin = CI_LL, ymax = CI_UL, fill = as.factor(model)), alpha = 0.1) +
  geom_line(aes(linetype = "estimated", col = as.factor(model))) + 
  geom_line(aes(y = Rt_real, linetype = "real")) + 
  facet_grid(cols = vars(dept_id), rows = vars(model2)) + 
  scale_color_viridis_d() + 
  scale_fill_viridis_d() +
  labs(linetype = "", col = "Rt assignment day", fill = "Rt assignment day") +
  theme_bw()
  
  
#### NPI lag ####
lag <- 0:9

cl <- makeCluster(8)
registerDoParallel(cl)

# Simulx
list_Rt_reg_lag_Simulx <- list()
list_Rt_reg_lag_fits_Simulx <- list()
for(i in 1:length(interval)){
  Rt_comp_res <- Rt_calc_fun(Inc_df = dataset1_Simulx, id_col_name = "dept_id", time_col_name = "day", 
                             Inc_col_name = "IncI_unscaled", model_name = lag[i], 
                             Rt_ref_df = true_Rt_df_Simulx, tstart = new_tstart[[7]], tend = new_tend[[7]], 
                             meansi = 10.1, stdsi = 8.75, Rt_prior = 1)
  
  reg <- Rt_reg_only_fun(data_for_est = Rt_comp_res$Rt_comp %>%
                           mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
                                  BG1 = ifelse(day > 70, 1, 0)), 
                         model_name = lag[i], lag_NPIs = TRUE, lag_days = lag[i], fits = TRUE)
  
  list_Rt_reg_lag_Simulx[[i]] <- reg$coefs
  list_Rt_reg_lag_fits_Simulx[[i]] <- reg$fits
  
}

reg_res_lag_Simulx_df <- do.call("rbind.data.frame", list_Rt_reg_lag_Simulx)
reg_res_lag_Simulx_fits_df <- do.call("rbind.data.frame", list_Rt_reg_lag_fits_Simulx)


list_Rt_reg_lag_Simulx_h <- list()
list_Rt_reg_lag_fits_Simulx_h <- list()
for(i in 1:length(interval)){
  Rt_comp_res <- Rt_calc_fun(Inc_df = dataset1_Simulx, id_col_name = "dept_id", time_col_name = "day", 
                             Inc_col_name = "IncH_unscaled", model_name = lag[i], 
                             Rt_ref_df = true_Rt_df_Simulx, tstart = new_tstart[[7]], tend = new_tend[[7]], 
                             meansi = 10.1, stdsi = 8.75, Rt_prior = 1)
  
  reg <- Rt_reg_only_fun(data_for_est = Rt_comp_res$Rt_comp %>%
                           mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
                                  BG1 = ifelse(day > 70, 1, 0)), 
                         model_name = lag[i], lag_NPIs = TRUE, lag_days = lag[i], fits = TRUE)
  
  list_Rt_reg_lag_Simulx_h[[i]] <- reg$coefs
  list_Rt_reg_lag_fits_Simulx_h[[i]] <- reg$fits
  
}

reg_res_lag_Simulx_h_df <- do.call("rbind.data.frame", list_Rt_reg_lag_Simulx_h)
reg_res_lag_Simulx_h_fits_df <- do.call("rbind.data.frame", list_Rt_reg_lag_fits_Simulx_h)


# ABM rm
list_Rt_reg_lag_ABM_rm <- list()
list_Rt_reg_lag_fits_ABM_rm <- list()
for(i in 1:length(interval)){
  Rt_comp_res <- Rt_calc_fun(Inc_df = data_ABM_rm_cov, id_col_name = "dept_id", time_col_name = "day", 
                             Inc_col_name = "IncI", model_name = lag[i], 
                             Rt_ref_df = true_Rt_df_ABM_rm, tstart = new_tstart[[7]], tend = new_tend[[7]], 
                             meansi = 8.2, stdsi = 5, Rt_prior = 2)

  reg <- Rt_reg_only_fun(data_for_est = Rt_comp_res$Rt_comp %>%
                           mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
                                  BG1 = ifelse(day > 70, 1, 0)), 
                         model_name = lag[i], lag_NPIs = TRUE, lag_days = lag[i], fits = TRUE)
  
  list_Rt_reg_lag_ABM_rm[[i]] <- reg$coefs
  list_Rt_reg_lag_fits_ABM_rm[[i]] <- reg$fits
}

reg_res_lag_ABM_rm_df <- do.call("rbind.data.frame", list_Rt_reg_lag_ABM_rm)
reg_res_lag_fits_ABM_rm_df <- do.call("rbind.data.frame", list_Rt_reg_lag_fits_ABM_rm)


# ABM hybrid
list_Rt_reg_lag_ABM_hybrid <- list()
list_Rt_reg_lag_fits_ABM_hybrid <- list()
for(i in 1:length(interval)){
  Rt_comp_res <- Rt_calc_fun(Inc_df = data_ABM_hybrid_cov, id_col_name = "dept_id", time_col_name = "day", 
                             Inc_col_name = "IncI", model_name = lag[i], 
                             Rt_ref_df = true_Rt_df_ABM_hybrid, tstart = new_tstart[[7]], tend = new_tend[[7]], 
                             meansi = 8.2, stdsi = 5, Rt_prior = 2)
  
  
  reg <- Rt_reg_only_fun(data_for_est = Rt_comp_res$Rt_comp %>%
                           mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
                                  BG1 = ifelse(day > 70, 1, 0)), 
                         model_name = lag[i], lag_NPIs = TRUE, lag_days = lag[i], fits = TRUE)
  
  list_Rt_reg_lag_ABM_hybrid[[i]] <- reg$coefs
  list_Rt_reg_lag_fits_ABM_hybrid[[i]] <- reg$fits
}

reg_res_lag_ABM_hybrid_df <- do.call("rbind.data.frame", list_Rt_reg_lag_ABM_hybrid)
reg_res_lag_fits_ABM_hybrid_df <- do.call("rbind.data.frame", list_Rt_reg_lag_fits_ABM_hybrid)


stopCluster(cl)

reg_comp_df_lag <- reg_res_lag_Simulx_df %>%
  mutate(model2 = "Simulx") %>%
  bind_rows(reg_res_lag_ABM_hybrid_df %>% mutate(model2 = "ABM hybrid")) %>%
  bind_rows(reg_res_lag_ABM_rm_df %>% mutate(model2 = "ABM rm")) %>%
  mutate(true_value = ifelse(parameter == "NPI 1", -1.45, -0.5))


ggplot(reg_comp_df_lag, aes(x = model, y = value, col = model2)) + 
  geom_pointrange(aes(ymin = CI_LL, ymax = CI_UL), position = position_dodge(width = 0.3)) + 
  geom_line(aes(y = true_value), linetype = "dashed", col = "darkred") +
  scale_x_continuous(breaks = 1:7) + 
  facet_wrap(~parameter) + 
  labs(y = "coefficient value", x = "lag", col = "model") + 
  scale_color_brewer(palette = "Dark2") +
  theme_bw()




ggplot(reg_res_lag_Simulx_h_fits_df %>% filter(dept_id %in% c(1, 4, 9, 13)), 
       aes(x = day, col = as.factor(lag))) + 
  geom_line(aes(y = Rt, linetype = "EpiEstim"), col = "black", linewidth = 1) +
  geom_line(aes(y = Rt_real, linetype = "Real"), col = "black", linewidth = 1) +
  geom_line(aes(y = Rt_fitted, linetype = "Reg fit"), linewidth = 0.8) +
  facet_wrap(~dept_id) + 
  labs(col = "NPI lag [days]", linetype = "Rt") + 
  scale_color_brewer(palette = "Dark2") +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) +
  theme_bw()




