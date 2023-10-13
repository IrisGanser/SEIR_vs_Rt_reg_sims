library(tidyverse)
library(lme4)
library(EpiEstim)
library(colorspace)
library(parallel)
library(foreach)
library(doParallel)


setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/sim_2params_regs")
dir <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_deSolve_data_creation_2params"
source("~/PhD/COVID_France/Dropbox_iris_covid/departement/Données_SPF/Data/data_functions.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/useful_functions.R")


popsize_df <- read.csv(paste0(dir, "/popsize_df.csv")) %>%
  mutate(dept_id = ifelse(dept_id > 20, dept_id-1, dept_id))

reg_res_list_I_2params_all <- vector(mode = "list")
reg_res_list_I_2params_random <- vector(mode = "list")
reg_res_list_I_2params_threshold <- vector(mode = "list")

reg_res_list_H_2params_all <- vector(mode = "list")
reg_res_list_H_2params_random <- vector(mode = "list")
reg_res_list_H_2params_threshold <- vector(mode = "list")

cl <- makeCluster(10)
registerDoParallel(cl)

for(j in 1:100){

    reg_data_all <- read.table(paste0(dir, "/data_sim_SEIRAHD_2params_all_ME", j, ".txt"), 
                           header = TRUE, sep = ",") %>%
      pivot_wider(names_from = obs_id, values_from = obs, names_prefix = "obs_") %>%
      rename(IncI = obs_3, IncH = obs_1)
    reg_data_threshold <- read.table(paste0(dir, "/data_sim_SEIRAHD_2params_threshold_ME", j, ".txt"), 
                               header = TRUE, sep = ",") %>%
      pivot_wider(names_from = obs_id, values_from = obs, names_prefix = "obs_") %>%
      rename(IncI = obs_3, IncH = obs_1)
    reg_data_random <- read.table(paste0(dir, "/data_sim_SEIRAHD_2params_random_ME", j, ".txt"), 
                               header = TRUE, sep = ",") %>%
      pivot_wider(names_from = obs_id, values_from = obs, names_prefix = "obs_") %>%
      rename(IncI = obs_3, IncH = obs_1)
    
    # infections
    res_all_I <- EpiEstim_reg_fun(data_for_est = reg_data_all, Inc_name = "IncI", rep_num = j)
    res_threshold_I <- EpiEstim_reg_fun(data_for_est = reg_data_threshold, Inc_name = "IncI", rep_num = j)
    res_random_I <- EpiEstim_reg_fun(data_for_est = reg_data_random, Inc_name = "IncI", rep_num = j)
    
    # hospitalizations, lagged by 5 days (mean time from infection to hospitalization)
    res_all_H <- EpiEstim_reg_fun(data_for_est = reg_data_all, Inc_name = "IncH", rep_num = j, 
                                  lag_NPIs = TRUE, lag_days = 5)
    res_threshold_H <- EpiEstim_reg_fun(data_for_est = reg_data_threshold, Inc_name = "IncH", rep_num = j, 
                                        lag_NPIs = TRUE, lag_days = 5)
    res_random_H <- EpiEstim_reg_fun(data_for_est = reg_data_random, Inc_name = "IncH", rep_num = j, 
                                     lag_NPIs = TRUE, lag_days = 5)

  
  reg_res_list_I_2params_all[[j]] <- res_all_I
  reg_res_list_I_2params_random[[j]] <- res_random_I
  reg_res_list_I_2params_threshold[[j]] <- res_threshold_I
  
  reg_res_list_H_2params_all[[j]] <- res_all_H
  reg_res_list_H_2params_random[[j]] <- res_random_H
  reg_res_list_H_2params_threshold[[j]] <- res_threshold_H
}

stopCluster(cl)

reg_res_I_2params_all_df <- do.call("rbind.data.frame", reg_res_list_I_2params_all)
reg_res_I_2params_threshold_df <- do.call("rbind.data.frame", reg_res_list_I_2params_threshold)
reg_res_I_2params_random_df <- do.call("rbind.data.frame", reg_res_list_I_2params_random)

reg_res_H_2params_all_df <- do.call("rbind.data.frame", reg_res_list_H_2params_all)
reg_res_H_2params_threshold_df <- do.call("rbind.data.frame", reg_res_list_H_2params_threshold)
reg_res_H_2params_random_df <- do.call("rbind.data.frame", reg_res_list_H_2params_random)

save(reg_res_I_2params_all_df, file = "reg_res_I_2params_all_df.RData")
save(reg_res_I_2params_random_df, file = "reg_res_I_2params_random_df.RData")
save(reg_res_I_2params_threshold_df, file = "reg_res_I_2params_threshold_df.RData")

save(reg_res_H_2params_all_df, file = "reg_res_H_2params_all_df.RData")
save(reg_res_H_2params_random_df, file = "reg_res_H_2params_random_df.RData")
save(reg_res_H_2params_threshold_df, file = "reg_res_H_2params_threshold_df.RData")



### regs ABM ####
dir2 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once2"
reg_res_list_I_2params_all_ABM_rm <- vector(mode = "list")

cl <- makeCluster(10)
registerDoParallel(cl)

for(j in 1:100){
  
  reg_data_all <- read.table(paste0(dir2, "/data_covasim_rm_", j, ".txt"), 
                             header = TRUE, sep = ",")
  
  res_all_I <- EpiEstim_reg_fun(data_for_est = reg_data_all, Inc_name = "IncI", rep_num = j)
  
  
  reg_res_list_I_2params_all_ABM_rm[[j]] <- res_all_I
  
}

stopCluster(cl)

reg_res_I_2params_all_ABM_rm_df <- do.call("rbind.data.frame", reg_res_list_I_2params_all_ABM_rm)
save(reg_res_I_2params_all_ABM_rm_df, file = "reg_res_I_2params_all_ABM_rm_df.RData")




#### new ABM models ####
dir4 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/boot_sim_2params_ABM5"
reg_res_list_I_2params_all_ABM_rm5 <- vector(mode = "list")
reg_res_list_I_2params_all_ABM_hybrid5 <- vector(mode = "list")

cl <- makeCluster(10)
registerDoParallel(cl)

# random mixing
for(j in 1:100){
  
  reg_data_all <- read.table(paste0(dir4, "/data_covasim_rm5_", j, ".txt"), 
                             header = TRUE, sep = ",")
  
  res_all_I <- EpiEstim_reg_fun(data_for_est = reg_data_all, Inc_name = "IncI", rep_num = j, 
                                meansi = 7.68, stdsi = 4.17)
  
  
  reg_res_list_I_2params_all_ABM_rm5[[j]] <- res_all_I
  
}


reg_res_I_2params_all_ABM_rm5_df <- do.call("rbind.data.frame", reg_res_list_I_2params_all_ABM_rm5)
save(reg_res_I_2params_all_ABM_rm5_df, file = "reg_res_I_2params_all_ABM_rm5_df.RData")


# hybrid
for(j in 1:99){
  
  reg_data_all <- read.table(paste0(dir4, "/data_covasim_hybrid5_", j, ".txt"), 
                             header = TRUE, sep = ",")
  
  res_all_I <- EpiEstim_reg_fun(data_for_est = reg_data_all, Inc_name = "IncI", rep_num = j, 
                                meansi = 7.68, stdsi = 4.17)
  
  
  reg_res_list_I_2params_all_ABM_hybrid5[[j]] <- res_all_I
  
}

stopCluster(cl)

reg_res_I_2params_all_ABM_hybrid5_df <- do.call("rbind.data.frame", reg_res_list_I_2params_all_ABM_hybrid5)
save(reg_res_I_2params_all_ABM_hybrid5_df, file = "reg_res_I_2params_all_ABM_hybrid5_df.RData")


# ABM 6 models
dir6 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/boot_sim_2params_ABM6"
reg_res_list_I_2params_all_ABM_rm6 <- vector(mode = "list")
reg_res_list_I_2params_all_ABM_hybrid6 <- vector(mode = "list")

cl <- makeCluster(10)
registerDoParallel(cl)

# random mixing
for(j in 1:100){
  
  reg_data_all <- read.table(paste0(dir6, "/data_SEIR_covasim_rm6_", j, ".txt"), 
                             header = TRUE, sep = ",")
  
  res_all_I <- EpiEstim_reg_fun(data_for_est = reg_data_all, Inc_name = "IncI", rep_num = j, 
                                meansi = 7.155, stdsi = 3.877)
  
  
  reg_res_list_I_2params_all_ABM_rm6[[j]] <- res_all_I
  
}


reg_res_I_2params_all_ABM_rm6_df <- do.call("rbind.data.frame", reg_res_list_I_2params_all_ABM_rm6)
save(reg_res_I_2params_all_ABM_rm6_df, file = "reg_res_I_2params_all_ABM_rm6_df.RData")


# hybrid
for(j in 1:100){
  
  reg_data_all <- read.table(paste0(dir6, "/data_SEIR_covasim_hybrid6_", j, ".txt"), 
                             header = TRUE, sep = ",")
  
  res_all_I <- EpiEstim_reg_fun(data_for_est = reg_data_all, Inc_name = "IncI", rep_num = j, 
                                meansi = 7.155, stdsi = 3.877)
  
  
  reg_res_list_I_2params_all_ABM_hybrid6[[j]] <- res_all_I
  
}

stopCluster(cl)

reg_res_I_2params_all_ABM_hybrid6_df <- do.call("rbind.data.frame", reg_res_list_I_2params_all_ABM_hybrid6)
save(reg_res_I_2params_all_ABM_hybrid6_df, file = "reg_res_I_2params_all_ABM_hybrid6_df.RData")



#### Simulx data creation
dir3 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"

reg_res_list_I_2params_all_Simulx <- vector(mode = "list")
reg_res_list_H_2params_all_Simulx <- vector(mode = "list")

cl <- makeCluster(10)
registerDoParallel(cl)

for(j in 1:100){
  
  reg_data_all <- read.table(paste0(dir3, "/data_sim_SEIRAHD_Simulx_2params_new2_ME", j, ".txt"), 
                             header = TRUE, sep = ",") %>%
    pivot_wider(names_from = obs_id, values_from = obs, names_prefix = "obs_") %>%
    rename(IncI = obs_3, IncH = obs_1)
  
  
  # infections
  res_all_I <- EpiEstim_reg_fun(data_for_est = reg_data_all, Inc_name = "IncI", rep_num = j, 
                                meansi = 10.1, stdsi = 8.75)
  
  # hospitalizations, lagged by 5 days (mean time from infection to hospitalization)
  res_all_H <- EpiEstim_reg_fun(data_for_est = reg_data_all, Inc_name = "IncH", rep_num = j, 
                                lag_NPIs = TRUE, lag_days = 5, 
                                meansi = 10.1, stdsi = 8.75)

  reg_res_list_I_2params_all_Simulx[[j]] <- res_all_I
  reg_res_list_H_2params_all_Simulx[[j]] <- res_all_H
  
}

stopCluster(cl)

reg_res_I_2params_all_Simulx_df <- do.call("rbind.data.frame", reg_res_list_I_2params_all_Simulx)
reg_res_H_2params_all_Simulx_df <- do.call("rbind.data.frame", reg_res_list_H_2params_all_Simulx)

save(reg_res_I_2params_all_Simulx_df, file = "reg_res_I_2params_all_Simulx_df.RData")
save(reg_res_H_2params_all_Simulx_df, file = "reg_res_H_2params_all_Simulx_df.RData")
