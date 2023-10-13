library(tidyverse)
library(lme4)
library(EpiEstim)
library(colorspace)
library(parallel)
library(foreach)
library(doParallel)


setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/sim_2params_regs")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")

#### Simulx SEIRAHD 2 models ####
dir2 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"


reg_res_list_I_Simulx_boot <- vector(mode = "list")
reg_res_list_H_Simulx_boot <- vector(mode = "list")

seeds <- seq(101, 600)

cl <- makeCluster(10)
registerDoParallel(cl)

# initialize final result lists
boot_res_3_levels <- list()
boot_res_3_levels_H <- list()


for(i in 1:100){ # i is the simulation dataset
 
  # initialize bootstrap lists
  reg_res_list_I_Simulx_boot <- list()
  reg_res_list_H_Simulx_boot <- list()
  
  
  # data for estimation
  reg_data_all <- read.table(paste0(dir2, "/data_sim_SEIRAHD_Simulx_2params_new2_ME", i, ".txt"), 
                             header = TRUE, sep = ",") %>%
    pivot_wider(names_from = obs_id, values_from = obs, names_prefix = "obs_") %>%
    rename(IncI = obs_3, IncH = obs_1) 
  
  
  for(j in 1:100){ # j is the bootstrap replicate
    
    # bootstrapping departments 
    set.seed(seeds[j])
    random_depts <- sample(1:94, 94, replace = TRUE)
    
    # Rt estimation
    Rt_list <- foreach(x = 1:94, .packages = c("tidyverse", "EpiEstim")) %dopar% {
      Inc_series <- reg_data_all %>% 
        filter(dept_id == random_depts[x]) 
      
      Inc <- as.numeric(na.omit(Inc_series$IncI*Inc_series$popsize/10^4))
      
      Rt_estim <- estimate_R(Inc,
                             method = "parametric_si",
                             config = make_config(list(
                               mean_si = 10.1,
                               std_si = 8.75)))$R
      Rt_estim <- Rt_estim %>%
        mutate(t = (t_end+t_start)/2,
               id = x, .before = t_start) %>%
        dplyr::select(t, id, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`) %>%
        rename(Rt = `Mean(R)`, CI_LL = `Quantile.0.025(R)`, CI_UL = `Quantile.0.975(R)`)
    }
    
    Rt_list_H <- foreach(x = 1:94, .packages = c("tidyverse", "EpiEstim")) %dopar% {
      Inc_series <- reg_data_all %>% 
        filter(dept_id == random_depts[x]) 
      
      Inc <- as.numeric(na.omit(Inc_series$IncH*Inc_series$popsize/10^4))
      
      Rt_estim <- estimate_R(Inc,
                             method = "parametric_si",
                             config = make_config(list(
                               mean_si = 10.1,
                               std_si = 8.75)))$R
      Rt_estim <- Rt_estim %>%
        mutate(t = (t_end+t_start)/2,
               id = x, .before = t_start) %>%
        dplyr::select(t, id, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`) %>%
        rename(Rt = `Mean(R)`, CI_LL = `Quantile.0.025(R)`, CI_UL = `Quantile.0.975(R)`)
    }
    
    
    Rt_df <- do.call("rbind.data.frame", Rt_list) %>%
      rename(dept_id = id, day = t) %>%
      mutate(se = (Rt - CI_LL)/1.96)  %>%
      mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
             BG1 = ifelse(day > 70, 1, 0))
    
    Rt_df_H <- do.call("rbind.data.frame", Rt_list_H) %>%
      rename(dept_id = id, day = t) %>%
      mutate(se = (Rt - CI_LL)/1.96) %>%
      mutate(lockdown1 = ifelse(between(day, 16, 70), 1, 0), 
             BG1 = ifelse(day > 70, 1, 0))
    
    
    
    # include uncertainty in Rt estimation by repeatedly sampling from its distribution
    reg_boot <- foreach(k = 1:500, .packages = c("tidyverse", "lme4")) %dopar% { # k is the Rt sampling
      # randomly sample one value for Rt from its distribution
      set.seed(seeds[k])
      Rt_df %<>% 
        filter(!is.na(Rt)) %>%
        rowwise() %>%
        mutate(Rt_sampled = rnorm(1, mean = Rt, sd = se)) %>%
        ungroup()
      
      # regression
      reg_res <- lmer(log(Rt_sampled) ~ lockdown1 + BG1 + (1|dept_id), data = Rt_df)
      
      
      # get reg coefficients
      coefs_Rt_reg <- coefficients(reg_res)$dept_id %>%
        dplyr::select(-1) %>%
        unique() %>%
        pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
        mutate(boot_rep = k)
      
    }
    
    reg_boot_H <- foreach(k = 1:500, .packages = c("tidyverse", "lme4")) %dopar% {
      # randomly sample one value for Rt from its distribution
      set.seed(seeds[k])
      Rt_df_H %<>% 
        filter(!is.na(Rt)) %>%
        rowwise() %>%
        mutate(Rt_sampled = rnorm(1, mean = Rt, sd = se)) %>%
        ungroup()
      
      # regression
      reg_res_H <- lmer(log(Rt_sampled) ~ lockdown1 + BG1 + (1|dept_id), data = Rt_df_H)
      
      # get reg coefficients
      coefs_Rt_reg_H <- coefficients(reg_res_H)$dept_id %>%
        dplyr::select(-1) %>%
        unique() %>%
        pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
        mutate(boot_rep = k)
    }
    
    
    # Regression results from one bootstrap sampling
    
    # point estimate regression (Rt without sampling)
    reg_point_est <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = Rt_df)
    reg_point_est_H <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = Rt_df_H)
    
    
    point_est <- coefficients(reg_point_est)$dept_id %>%
      dplyr::select(-1) %>%
      unique() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") 
    
    point_est_H <- coefficients(reg_point_est_H)$dept_id %>%
      dplyr::select(-1) %>%
      unique() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") 
    
    
    # piece point estimate and CI together
    reg_boot_df <- do.call("rbind.data.frame", reg_boot) %>%
      group_by(parameter) %>%
      summarize(CI_LL = quantile(value, probs = 0.025), 
                CI_UL = quantile(value, probs = 0.975)) %>%
      left_join(point_est, by = "parameter") %>%
      mutate(rep = j, sim_rep = i)
    
    reg_boot_df_H <- do.call("rbind.data.frame", reg_boot_H) %>%
      group_by(parameter) %>%
      summarize(CI_LL = quantile(value, probs = 0.025), 
                CI_UL = quantile(value, probs = 0.975)) %>%
      left_join(point_est_H, by = "parameter") %>%
      mutate(rep = j, sim_rep = i)
    
    
    reg_res_list_I_Simulx_boot[[j]] <- reg_boot_df
    reg_res_list_H_Simulx_boot[[j]] <- reg_boot_df_H

  }
  
  # summarize results from 100 bootstrap runs into a data frame
  boot_Simulx_I_df <- do.call("rbind.data.frame", reg_res_list_I_Simulx_boot)
  boot_Simulx_H_df <- do.call("rbind.data.frame", reg_res_list_H_Simulx_boot)
  
  # save data frames in the final lists
  boot_res_3_levels[[i]] <- boot_Simulx_I_df
  boot_res_3_levels_H[[i]] <- boot_Simulx_H_df
}

df_boot_res_3_levels <- do.call("rbind.data.frame", boot_res_3_levels)
df_boot_res_3_levels_H <- do.call("rbind.data.frame", boot_res_3_levels_H)

save(df_boot_res_3_levels, file = "df_boot_res_3_levels.RData")
save(df_boot_res_3_levels_H, file = "df_boot_res_3_levels_H.RData")