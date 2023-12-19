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
library(lixoftConnectors)
initializeLixoftConnectors(software = "simulx")

setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/Rt_trajectories")
source("~/PhD/COVID_France/Dropbox_iris_covid/departement/Données_SPF/Data/data_functions.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")


dir1 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"
source(paste0(dir1, "/SEIRAHD_Simulx_function_2params.R"))
  


# load initial parameters
load(paste(dir1, "init_list.RData", sep = "/"))




# Simulx 2 ----------------------------------------------------------------


#### modify regression parameters ####
popsize_df <- read.csv(paste(dir1, "popsize_df.csv", sep = "/")) %>%
  mutate(dept_id = ifelse(dept_id < 20, dept_id, dept_id - 1))


# ld1 start
ld1_start <- c(20, 30, 40, 50, 60)

for(i in 1:length(ld1_start)){
  ld1_reg_df <- data.frame(id = rep(1:94, each = 151), 
                           time = rep(1:151, 94), 
                           initH = rep(init_list[[1]]$H, each = 151)) %>%
    mutate(lockdown1 = ifelse(between(time, ld1_start[i], 50 + ld1_start[i]), 1, 0),
           BG1 = ifelse(time > 50 + ld1_start[i], 1, 0), 
           popsize = 10000)
  
  write.csv(ld1_reg_df, paste0("ld1_reg_df_2params_late", i, ".csv"), row.names = F) 
}

# individual parameters
ind_params <- read.table(paste0(dir1, "/ind_params/ind_2params_1.txt"), sep = " ", header = TRUE) %>%
  dplyr::select(-c(beta_school, beta_w)) 

ind_params_lowb1 <- ind_params %>%
  mutate(b1 = ifelse(b1 > 0.4, b1 - 0.15, b1))


write.table(ind_params, file = "ind_params.txt", sep = " ", row.names = F)
write.table(ind_params_lowb1, file = "ind_params_lowb1.txt", sep = " ", row.names = F)


ld1_strength <- seq(0.5, 2, 0.1)
for(i in 1:length(ld1_strength)){
  ind_params_x <- ind_params %>%
    mutate(beta_ld1 = ld1_strength[i])
  
  write.table(ind_params_x, paste0("ind_params", ld1_strength[i], ".txt"), 
              sep = " ", row.names = F)
}

#### Simulx data simulation ####
project.file <- paste(dir1, "sim_SEIRAHD_2params_init_est.mlxtran", sep = "/")
importMonolixProject(project.file)

# ld1 start
sim_res_ld1_start <- list()

for(j in 1:length(ld1_start)){
  
  sim_res <- sim_SEIRAHD_Simulx_ME_2params(path_to_ind_params = paste0(getwd(), "/ind_params_lowb1.txt"), 
                                           regressor_df_path = paste0("ld1_reg_df_2params_late", j, ".csv"), 
                                           tmax = 151)
  
  sim_res_ld1_start[[j]] <- sim_res
  
  data_for_sim <- read.csv(paste0("ld1_reg_df_2params_late", j, ".csv"))
  
  monolix_SEIRAHD <- data_for_sim %>%
    left_join(., sim_res %>% 
                select(id, time, IncI_ME, IncH_ME, PrevH_ME, IncD_ME), 
              by = c("id", "time")) %>%
    rename(dept_id = id, day = time) %>% 
    group_by(dept_id) %>%
    mutate(initH = PrevH_ME[day == 1], .before = lockdown1) %>%
    ungroup() %>%
    pivot_longer(c(IncI_ME, IncH_ME, PrevH_ME, IncD_ME), 
                 names_to = "obs_id", values_to = "obs") %>%
    mutate(obs_id = case_when(obs_id == "IncH_ME" ~ 1, 
                              obs_id == "PrevH_ME" ~ 2,
                              obs_id == "IncI_ME" ~ 3, 
                              obs_id == "IncD_ME" ~ 4)) %>%
    relocate(obs, .after = day) %>%
    relocate(obs_id, .after = obs) %>%
    mutate(obs = ifelse(obs < 0, 0, obs)) # correct data in case some observations are < 0
  
  write.table(monolix_SEIRAHD, file = paste0("data_sim_SEIRAHD_2params_ld1_start", ld1_start[j], ".txt"), 
              row.names = FALSE, sep = ",")
  
  monolix_SEIR <- monolix_SEIRAHD %>%
    filter(obs_id == 3) %>%
    rename(IncI = obs) %>%
    select(-c(obs_id, initH))
  
  write.table(monolix_SEIR, file = paste0("data_sim_SEIR_Simulx_2params_ld1_start", ld1_start[j], ".txt"), 
              sep = ",", row.names = FALSE)
}

save(sim_res_ld1_start, file = "sim_res_ld1_start.RData")


for(j in 1:length(ld1_start)){
  data <- read.table(paste0("data_sim_SEIR_Simulx_2params_ld1_start", ld1_start[j], ".txt"), 
                     sep = ",", header = TRUE)
  
  plot <- ggplot(data, aes(x = day, y = IncI, group = dept_id)) + 
    geom_line() +
    labs(title = paste("Sim Ld1 start", ld1_start[j]))

  
  print(plot)
}


# ld1 strength
sim_res_ld1_strength <- list()

for(j in 1:length(ld1_strength)){
  
  sim_res <- sim_SEIRAHD_Simulx_ME_2params(path_to_ind_params = paste0(getwd(), "/ind_params", ld1_strength[j], ".txt"), 
                                           regressor_df_path = "ld1_reg_df_2params_late1.csv", 
                                           tmax = 151)
  
  sim_res_ld1_strength[[j]] <- sim_res
  
  
  data_for_sim <- read.csv("ld1_reg_df_2params_late1.csv")
  
  monolix_SEIRAHD <- data_for_sim %>%
    left_join(., sim_res %>% 
                select(id, time, IncI_ME, IncH_ME, PrevH_ME, IncD_ME), 
              by = c("id", "time")) %>%
    rename(dept_id = id, day = time) %>% 
    group_by(dept_id) %>%
    mutate(initH = PrevH_ME[day == 1], .before = lockdown1) %>%
    ungroup() %>%
    pivot_longer(c(IncI_ME, IncH_ME, PrevH_ME, IncD_ME), 
                 names_to = "obs_id", values_to = "obs") %>%
    mutate(obs_id = case_when(obs_id == "IncH_ME" ~ 1, 
                              obs_id == "PrevH_ME" ~ 2,
                              obs_id == "IncI_ME" ~ 3, 
                              obs_id == "IncD_ME" ~ 4)) %>%
    relocate(obs, .after = day) %>%
    relocate(obs_id, .after = obs) %>%
    mutate(obs = ifelse(obs < 0, 0, obs)) # correct data in case some observations are < 0
  
  write.table(monolix_SEIRAHD, file = paste0("data_sim_SEIRAHD_2params_ld1_strength", ld1_strength[j], ".txt"), 
              row.names = FALSE, sep = ",")
  
  monolix_SEIR <- monolix_SEIRAHD %>%
    filter(obs_id == 3) %>%
    rename(IncI = obs) %>%
    select(-c(obs_id, initH))
  
  write.table(monolix_SEIR, file = paste0("data_sim_SEIR_Simulx_2params_ld1_strength", ld1_strength[j], ".txt"), 
              sep = ",", row.names = FALSE)
}

save(sim_res_ld1_strength, file = "sim_res_ld1_strength.RData")


for(j in 1:length(ld1_strength)){
  data <- read.table(paste0("data_sim_SEIR_Simulx_2params_ld1_strength", ld1_strength[j], ".txt"), 
                     sep = ",", header = TRUE)
  
  plot <- ggplot(data, aes(x = day, y = IncI, group = dept_id)) + 
    geom_line() +
    labs(title = paste("Sim Ld1 strength", ld1_strength[j]))
  
  print(plot)
}


# higher basic transmission
sim_res_ld1_start_high <- list()

for(j in 1:length(ld1_start)){
  
  sim_res <- sim_SEIRAHD_Simulx_ME_2params(path_to_ind_params = "ind_params.txt", 
                                           regressor_df_path = paste0("ld1_reg_df_2params_late", j, ".csv"), 
                                           tmax = 151)
  
  sim_res_ld1_start_high[[j]] <- sim_res
  
  data_for_sim <- read.csv(paste0("ld1_reg_df_2params_late", j, ".csv"))
  
  monolix_SEIRAHD <- data_for_sim %>%
    left_join(., sim_res %>% 
                select(id, time, IncI_ME, IncH_ME, PrevH_ME, IncD_ME), 
              by = c("id", "time")) %>%
    rename(dept_id = id, day = time) %>% 
    group_by(dept_id) %>%
    mutate(initH = PrevH_ME[day == 1], .before = lockdown1) %>%
    ungroup() %>%
    pivot_longer(c(IncI_ME, IncH_ME, PrevH_ME, IncD_ME), 
                 names_to = "obs_id", values_to = "obs") %>%
    mutate(obs_id = case_when(obs_id == "IncH_ME" ~ 1, 
                              obs_id == "PrevH_ME" ~ 2,
                              obs_id == "IncI_ME" ~ 3, 
                              obs_id == "IncD_ME" ~ 4)) %>%
    relocate(obs, .after = day) %>%
    relocate(obs_id, .after = obs) %>%
    mutate(obs = ifelse(obs < 0, 0, obs)) # correct data in case some observations are < 0
  
  
  write.table(monolix_SEIRAHD, file = paste0("data_sim_SEIRAHD_2params_ld1_start_high_", ld1_start[j], ".txt"), 
              row.names = FALSE, sep = ",")
  
  monolix_SEIR <- monolix_SEIRAHD %>%
    filter(obs_id == 3) %>%
    rename(IncI = obs) %>%
    select(-c(obs_id, initH))
  
  write.table(monolix_SEIR, file = paste0("data_sim_SEIR_Simulx_2params_ld1_start_high_", ld1_start[j], ".txt"), 
              sep = ",", row.names = FALSE)
}

save(sim_res_ld1_start_high, file = "sim_res_ld1_start_high.RData")


for(j in 1:length(ld1_start)){
  data <- read.table(paste0("data_sim_SEIR_Simulx_2params_ld1_start_high_", ld1_start[j], ".txt"), 
                     sep = ",", header = TRUE)
  
  plot <- ggplot(data, aes(x = day, y = IncI, group = dept_id)) + 
    geom_line() +
    labs(title = paste("Sim Ld1 start high", ld1_start[j]))
  
  
  print(plot)
}



# Simulx 3 ----------------------------------------------------------------

# individual parameters
ind_params3 <- read.table(paste0(dir1, "/ind_params/ind_2params_1.txt"), sep = " ", header = TRUE) %>%
  dplyr::select(-c(beta_school, beta_w)) %>%
  mutate(beta_BG1 = 0.8)

ind_params3_lowb1 <- ind_params3 %>%
  mutate(b1 = ifelse(b1 > 0.4, b1 - 0.15, b1))

write.table(ind_params3, file = "ind_params3.txt", sep = " ", row.names = F)
write.table(ind_params3_lowb1, file = "ind_params3_lowb1.txt", sep = " ", row.names = F)


ld1_strength <- seq(0.5, 2, 0.1)
for(i in 1:length(ld1_strength)){
  ind_params_x <- ind_params3 %>%
    mutate(beta_ld1 = ld1_strength[i])
  
  write.table(ind_params_x, paste0("ind_params3_", ld1_strength[i], ".txt"), 
              sep = " ", row.names = F)
}

#### Simulx data simulation ####
project.file <- paste(dir1, "sim_SEIRAHD_2params_init_est.mlxtran", sep = "/")
importMonolixProject(project.file)

# ld1 start
sim_res_ld1_start3 <- list()

for(j in 1:length(ld1_start)){
  
  sim_res <- sim_SEIRAHD_Simulx_ME_2params(path_to_ind_params = paste0(getwd(), "/ind_params3_lowb1.txt"), 
                                           regressor_df_path = paste0("ld1_reg_df_2params_late", j, ".csv"), 
                                           tmax = 151)
  
  sim_res_ld1_start3[[j]] <- sim_res
  
  data_for_sim <- read.csv(paste0("ld1_reg_df_2params_late", j, ".csv"))
  
  monolix_SEIRAHD <- data_for_sim %>%
    left_join(., sim_res %>% 
                select(id, time, IncI_ME, IncH_ME, PrevH_ME, IncD_ME), 
              by = c("id", "time")) %>%
    rename(dept_id = id, day = time) %>% 
    group_by(dept_id) %>%
    mutate(initH = PrevH_ME[day == 1], .before = lockdown1) %>%
    ungroup() %>%
    pivot_longer(c(IncI_ME, IncH_ME, PrevH_ME, IncD_ME), 
                 names_to = "obs_id", values_to = "obs") %>%
    mutate(obs_id = case_when(obs_id == "IncH_ME" ~ 1, 
                              obs_id == "PrevH_ME" ~ 2,
                              obs_id == "IncI_ME" ~ 3, 
                              obs_id == "IncD_ME" ~ 4)) %>%
    relocate(obs, .after = day) %>%
    relocate(obs_id, .after = obs) %>%
    mutate(obs = ifelse(obs < 0, 0, obs)) # correct data in case some observations are < 0
  
  write.table(monolix_SEIRAHD, file = paste0("data_sim_SEIRAHD_2params3_ld1_start", ld1_start[j], ".txt"), 
              row.names = FALSE, sep = ",")
  
  monolix_SEIR <- monolix_SEIRAHD %>%
    filter(obs_id == 3) %>%
    rename(IncI = obs) %>%
    select(-c(obs_id, initH))
  
  write.table(monolix_SEIR, file = paste0("data_sim_SEIR_Simulx_2params3_ld1_start", ld1_start[j], ".txt"), 
              sep = ",", row.names = FALSE)
}

save(sim_res_ld1_start3, file = "sim_res_params3_ld1_start.RData")


for(j in 1:length(ld1_start)){
  data <- read.table(paste0("data_sim_SEIR_Simulx_2params3_ld1_start", ld1_start[j], ".txt"), 
                     sep = ",", header = TRUE)
  
  plot <- ggplot(data, aes(x = day, y = IncI, group = dept_id)) + 
    geom_line() +
    labs(title = paste("Sim Ld1 start", ld1_start[j]))
  
  
  print(plot)
}


# ld1 strength
sim_res_ld1_strength3 <- list()

for(j in 1:length(ld1_strength)){
  
  sim_res <- sim_SEIRAHD_Simulx_ME_2params(path_to_ind_params = paste0(getwd(), "/ind_params3_", ld1_strength[j], ".txt"), 
                                           regressor_df_path = "ld1_reg_df_2params_late1.csv", 
                                           tmax = 151)
  
  sim_res_ld1_strength3[[j]] <- sim_res
  
  
  data_for_sim <- read.csv("ld1_reg_df_2params_late1.csv")
  
  monolix_SEIRAHD <- data_for_sim %>%
    left_join(., sim_res %>% 
                select(id, time, IncI_ME, IncH_ME, PrevH_ME, IncD_ME), 
              by = c("id", "time")) %>%
    rename(dept_id = id, day = time) %>% 
    group_by(dept_id) %>%
    mutate(initH = PrevH_ME[day == 1], .before = lockdown1) %>%
    ungroup() %>%
    pivot_longer(c(IncI_ME, IncH_ME, PrevH_ME, IncD_ME), 
                 names_to = "obs_id", values_to = "obs") %>%
    mutate(obs_id = case_when(obs_id == "IncH_ME" ~ 1, 
                              obs_id == "PrevH_ME" ~ 2,
                              obs_id == "IncI_ME" ~ 3, 
                              obs_id == "IncD_ME" ~ 4)) %>%
    relocate(obs, .after = day) %>%
    relocate(obs_id, .after = obs) %>%
    mutate(obs = ifelse(obs < 0, 0, obs)) # correct data in case some observations are < 0
  
  write.table(monolix_SEIRAHD, file = paste0("data_sim_SEIRAHD_2params3_ld1_strength", ld1_strength[j], ".txt"), 
              row.names = FALSE, sep = ",")
  
  monolix_SEIR <- monolix_SEIRAHD %>%
    filter(obs_id == 3) %>%
    rename(IncI = obs) %>%
    select(-c(obs_id, initH))
  
  write.table(monolix_SEIR, file = paste0("data_sim_SEIR_Simulx_2params3_ld1_strength", ld1_strength[j], ".txt"), 
              sep = ",", row.names = FALSE)
}

save(sim_res_ld1_strength3, file = "sim_res_params3_ld1_strength.RData")


for(j in 1:length(ld1_strength)){
  data <- read.table(paste0("data_sim_SEIR_Simulx_2params3_ld1_strength", ld1_strength[j], ".txt"), 
                     sep = ",", header = TRUE)
  
  plot <- ggplot(data, aes(x = day, y = IncI, group = dept_id)) + 
    geom_line() +
    labs(title = paste("Sim Ld1 strength", ld1_strength[j]))
  
  print(plot)
}


# higher basic transmission
sim_res_ld1_start_high3 <- list()

for(j in 1:length(ld1_start)){
  
  sim_res <- sim_SEIRAHD_Simulx_ME_2params(path_to_ind_params = "ind_params3.txt", 
                                           regressor_df_path = paste0("ld1_reg_df_2params_late", j, ".csv"), 
                                           tmax = 151)
  
  sim_res_ld1_start_high3[[j]] <- sim_res
  
  data_for_sim <- read.csv(paste0("ld1_reg_df_2params_late", j, ".csv"))
  
  monolix_SEIRAHD <- data_for_sim %>%
    left_join(., sim_res %>% 
                select(id, time, IncI_ME, IncH_ME, PrevH_ME, IncD_ME), 
              by = c("id", "time")) %>%
    rename(dept_id = id, day = time) %>% 
    group_by(dept_id) %>%
    mutate(initH = PrevH_ME[day == 1], .before = lockdown1) %>%
    ungroup() %>%
    pivot_longer(c(IncI_ME, IncH_ME, PrevH_ME, IncD_ME), 
                 names_to = "obs_id", values_to = "obs") %>%
    mutate(obs_id = case_when(obs_id == "IncH_ME" ~ 1, 
                              obs_id == "PrevH_ME" ~ 2,
                              obs_id == "IncI_ME" ~ 3, 
                              obs_id == "IncD_ME" ~ 4)) %>%
    relocate(obs, .after = day) %>%
    relocate(obs_id, .after = obs) %>%
    mutate(obs = ifelse(obs < 0, 0, obs)) # correct data in case some observations are < 0
  
  
  write.table(monolix_SEIRAHD, file = paste0("data_sim_SEIRAHD_2param3_ld1_start_high_", ld1_start[j], ".txt"), 
              row.names = FALSE, sep = ",")
  
  monolix_SEIR <- monolix_SEIRAHD %>%
    filter(obs_id == 3) %>%
    rename(IncI = obs) %>%
    select(-c(obs_id, initH))
  
  write.table(monolix_SEIR, file = paste0("data_sim_SEIR_Simulx_2params3_ld1_start_high_", ld1_start[j], ".txt"), 
              sep = ",", row.names = FALSE)
}

save(sim_res_ld1_start_high3, file = "sim_res_params3_ld1_start_high.RData")


for(j in 1:length(ld1_start)){
  data <- read.table(paste0("data_sim_SEIR_Simulx_2params3_ld1_start_high_", ld1_start[j], ".txt"), 
                     sep = ",", header = TRUE)
  
  plot <- ggplot(data, aes(x = day, y = IncI, group = dept_id)) + 
    geom_line() +
    labs(title = paste("Sim Ld1 start high", ld1_start[j]))
  
  
  print(plot)
}
