EpiEstim_reg_fun <- function(data_for_est, Inc_name, rep_num, lag_NPIs = FALSE, lag_days = 0, 
                             meansi = 7.5, stdsi = 5, meanprior = 2, stdprior = 4){

  
  Rt_list <- foreach(i = 1:94, .packages = c("tidyverse", "EpiEstim")) %dopar% {
    Inc_series <- data_for_est %>% 
      filter(dept_id == i) 
    
    Inc <- as.numeric(na.omit(Inc_series[[Inc_name]]*Inc_series$popsize/10^4))
    
    Rt_estim <- estimate_R(Inc,
                           method = "parametric_si", 
                           config = make_config(list(
                             mean_si = meansi,
                             std_si = stdsi, 
                             mean_prior = meanprior, 
                             std_prior = stdprior)))$R
    Rt_estim <- Rt_estim %>%
      mutate(t = (t_end+t_start)/2,
             id = i, .before = t_start) %>%
      dplyr::select(t, id, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`) %>%
      rename(Rt = `Mean(R)`, CI_LL = `Quantile.0.025(R)`, CI_UL = `Quantile.0.975(R)`)
  }
  
  Rt_df <- do.call("rbind.data.frame", Rt_list) %>%
    rename(dept_id = id, day = t)
  
  if(lag_NPIs){
    
    reg_data <- data_for_est %>%
      dplyr::select(dept_id, day, lockdown1, BG1) %>%
      group_by(dept_id) %>%
      mutate(lockdown1 = lag(lockdown1, n = lag_days, default = 0), 
             BG1 = lag(BG1, n = lag_days, default = 0)) %>%
      ungroup() %>%
      unique() %>%
      left_join(., Rt_df, by = c("dept_id", "day"))
    
    reg_res <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data)
    
  } else{
    
    reg_data <- data_for_est %>%
      dplyr::select(dept_id, day, lockdown1, BG1) %>%
      unique() %>%
      left_join(., Rt_df, by = c("dept_id", "day"))
    
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
                              labels = c("Lockdown 1", "Barrier gestures"))) %>%
    mutate(rep = rep_num)
  
  return(coefs_Rt_reg)
  
}


popparam_cleaning <- function(df){
  df_clean <- df %>%
    mutate(parameter = str_remove(parameter, "_pop"), 
           value = 0 - value) %>% # invert the sign of the parameter because it's estimated different in Monolix
    filter(grepl("beta", parameter) & !grepl("omega", parameter)) 
  
  return(df_clean)
}

bootstrap_cleaning <- function(df, boot_rep = i, sim_rep = j){
  df_clean <- df %>%
    mutate(parameter = str_remove(parameter, "_pop"), 
          value = 0 - value) %>% # invert the sign of the parameter because it's estimated different in Monolix
    filter(grepl("beta", parameter) & !grepl("omega", parameter)) %>%
    mutate(boot_rep = boot_rep, sim_rep = sim_rep)
  
  return(df_clean)
}

bootstrap_CI_calc <- function(boot_list){
  boot_df <- do.call("rbind.data.frame", boot_list) %>%
    group_by(parameter, sim_rep) %>%
    summarize(mean_est1 = mean(value), 
              sd_est = sd(value), 
              CI_LL2 = quantile(value, probs = 0.025), 
              CI_UL2 = quantile(value, probs = 0.975)) %>%
    mutate(CI_LL1 = mean_est1 - 1.96*sd_est, 
           CI_UL1 = mean_est1 + 1.96*sd_est)
  
  return(boot_df)
}


bootstrap_summary <- function(bootstrap_list, true_val_NPI1 = -1.45, true_val_NPI2 = -0.5){
  boot_df <- do.call("rbind.data.frame", bootstrap_list) %>%
    ungroup() %>% 
    mutate(parameter = factor(parameter, 
                              levels = c("beta_ld1", "beta_BG1"),
                              labels = c("NPI 1", "NPI 2")), 
           true_value = ifelse(parameter == "NPI 1", true_val_NPI1, true_val_NPI2), 
           unique_sims = n()/2, 
           CI_covers = ifelse(between(true_value, CI_LL2, CI_UL2), 1, 0), 
           bias = true_value - mean_est2, 
           rel_bias = abs(true_value - mean_est2)/abs(true_value)*100) %>%
    group_by(parameter) %>%
    mutate(perc_CI_covers = sum(CI_covers)/unique_sims*100, 
           mean_bias = mean(bias, na.rm = TRUE), 
           mean_rel_bias = mean(rel_bias, na.rm = TRUE))
  
  return(boot_df)
}


reg_summary <- function(reg_df, true_val_NPI1 = -1.45, true_val_NPI2 = -0.5){
  reg_df %<>% 
    mutate(parameter = factor(parameter, 
                              levels = c("Lockdown 1", "Barrier gestures"), 
                              labels = c("NPI 1", "NPI 2")), 
           true_value = ifelse(parameter == "NPI 1", true_val_NPI1, true_val_NPI2), 
           unique_sims = length(unique(rep)), 
           CI_covers = ifelse(between(true_value, CI_LL, CI_UL), 1, 0), 
           bias = true_value - value, 
           rel_bias = abs(true_value - value)/abs(true_value)*100) %>%
    group_by(parameter, model) %>%
    mutate(perc_CI_covers = sum(CI_covers)/unique_sims*100, 
              mean_bias = mean(bias), 
              mean_rel_bias = mean(rel_bias))
    
  return(reg_df)
}


SEIR_summary_2params <- function(SEIR_list){
  df <- do.call("rbind.data.frame", SEIR_list) %>%
    mutate(parameter = str_remove(parameter, "_pop"))
  
  df_boot <- df %>%
    group_by(parameter, model) %>%
    mutate(value = 0-value) %>%
    summarize(median_est = median(value),
              CI_LL = quantile(value, probs = 0.025),
              CI_UL = quantile(value, probs = 0.975)) %>%
    mutate(parameter = str_remove(parameter, "_pop")) %>%
    filter(grepl("beta", parameter)) %>%
    filter(!grepl("omega", parameter)) %>%
    mutate(parameter = str_remove(parameter, "beta_")) %>%
    mutate(parameter = factor(parameter, 
                              levels = c("ld1", "BG1"), 
                              labels = c("Lockdown 1", "Barrier gestures"))) %>%
    rename(value = median_est)
}


SEIR_summary_2params2 <- function(SEIR_list){
  df <- do.call("rbind.data.frame", SEIR_list) %>%
    mutate(parameter = str_remove(parameter, "_pop")) %>%
    filter(grepl("beta", parameter) & !grepl("omega", parameter)) %>%
    mutate(parameter = str_remove(parameter, "beta_")) %>%
    mutate(parameter = factor(parameter, 
                              levels = c("ld1", "BG1"), 
                              labels = c("Lockdown 1", "Barrier gestures")), 
           value = 0-value)
}


calc_Rt_SEIR <- function(b1, gamma, S, N){
  R_eff <- b1*S/(gamma*N)
  return(R_eff)
}


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


EpiEstim_only_fun <- function(data_for_est, Inc_name, meansi = 7.5, stdsi = 5, 
                              meanprior = 2, stdprior = 4){
  
  
  Rt_list <- foreach(i = 1:94, .packages = c("tidyverse", "EpiEstim")) %dopar% {
    Inc_series <- data_for_est %>% 
      filter(dept_id == i) 
    
    Inc <- as.numeric(na.omit(Inc_series[[Inc_name]]*Inc_series$popsize/10^4))
    
    Rt_estim <- estimate_R(Inc,
                           method = "parametric_si",
                           config = make_config(list(
                             mean_si = meansi,
                             std_si = stdsi, 
                             mean_prior = meanprior, 
                             str_prior = stdprior)))$R
    Rt_estim <- Rt_estim %>%
      mutate(t = (t_end+t_start)/2,
             id = i, .before = t_start) %>%
      dplyr::select(t, id, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`) %>%
      rename(Rt = `Mean(R)`, CI_LL = `Quantile.0.025(R)`, CI_UL = `Quantile.0.975(R)`)
  }
  
  Rt_df <- do.call("rbind.data.frame", Rt_list) %>%
    rename(dept_id = id, day = t)
  
  return(Rt_df)
}

