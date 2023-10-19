library(tidyverse)
library(deSolve)
library(lme4)
library(EpiEstim)
library(colorspace)
library(parallel)
library(foreach)
library(doParallel)
library(lixoftConnectors)
initializeLixoftConnectors(software = "simulx")

setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/Rt_trajectories")
source("~/PhD/COVID_France/Dropbox_iris_covid/departement/Données_SPF/Data/data_functions.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")


dir1 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"
source(paste0(dir1, "/SEIRAHD_Simulx_function_2params.R"))


# load initial parameters
load(paste(dir1, "init_list.RData", sep = "/"))

#### modify regression parameters ####
popsize_df <- read.csv(paste(dir1, "popsize_df.csv", sep = "/")) %>%
  mutate(dept_id = ifelse(dept_id < 20, dept_id, dept_id - 1))


ld1_reg_1week_slow <- data.frame(id = rep(1:94, each = 121), 
                                 time = rep(1:121, 94), 
                                 initH = rep(init_list[[1]]$H, each = 121)) %>%
  mutate(lockdown1 = rep(c(rep(0, 14), seq(0, 1, length.out = 8), rep(1, 70-23), 
                       seq(1, 0, length.out = 8), rep(0, 51-7)), 94),
         BG1 = rep(c(rep(0, 69), seq(0, 1, length.out = 8), rep(1, 121-69-8)), 94),
         popsize = 10000)

ld1_reg_2week_slow <- data.frame(id = rep(1:94, each = 121), 
                                 time = rep(1:121, 94), 
                                 initH = rep(init_list[[1]]$H, each = 121)) %>%
  mutate(lockdown1 = rep(c(rep(0, 14), seq(0, 1, length.out = 15), rep(1, 40), 
                           seq(1, 0, length.out = 15), rep(0, 51-14)), 94),
         BG1 = rep(c(rep(0, 69), seq(0, 1, length.out = 15), rep(1, 121-69-15)), 94),
         popsize = 10000)


ld1_reg_1week_slow_late <- data.frame(id = rep(1:94, each = 121), 
                                 time = rep(1:121, 94), 
                                 initH = rep(init_list[[1]]$H, each = 121)) %>%
  mutate(lockdown1 = rep(c(rep(0, 25), seq(0, 1, length.out = 8), rep(1, 50), 
                           seq(1, 0, length.out = 8), rep(0, 30)), 94),
         BG1 = rep(c(rep(0, 83), seq(0, 1, length.out = 8), rep(1, 121-83-8)), 94),
         popsize = 10000)

ld1_reg_2week_slow_late <- data.frame(id = rep(1:94, each = 121), 
                                 time = rep(1:121, 94), 
                                 initH = rep(init_list[[1]]$H, each = 121)) %>%
  mutate(lockdown1 = rep(c(rep(0, 25), seq(0, 1, length.out = 15), rep(1, 40), 
                           seq(1, 0, length.out = 15), rep(0, 26)), 94),
         BG1 = rep(c(rep(0, 79), seq(0, 1, length.out = 15), rep(1, 121-79-15)), 94),
         popsize = 10000)


write.csv(ld1_reg_1week_slow, "ld1_reg_1week_slow.csv", row.names = F) 
write.csv(ld1_reg_2week_slow, "ld1_reg_2week_slow.csv", row.names = F) 
write.csv(ld1_reg_1week_slow_late, "ld1_reg_1week_slow_late.csv", row.names = F) 
write.csv(ld1_reg_2week_slow_late, "ld1_reg_2week_slow_late.csv", row.names = F) 

list_reg_files_slow_NPI <- list(ld1_reg_1week_slow, ld1_reg_2week_slow, 
                                ld1_reg_1week_slow_late, ld1_reg_2week_slow_late)



#### Simulx data simulation ####
project.file <- paste(dir1, "sim_SEIRAHD_2params_init_est.mlxtran", sep = "/")
importMonolixProject(project.file)




sim_res_1week_slow <- sim_SEIRAHD_Simulx_ME_2params(path_to_ind_params = paste0(getwd(), "/ind_params.txt"), 
                                           regressor_df_path = "ld1_reg_1week_slow.csv")
sim_res_2week_slow <- sim_SEIRAHD_Simulx_ME_2params(path_to_ind_params = paste0(getwd(), "/ind_params.txt"), 
                                                    regressor_df_path = "ld1_reg_2week_slow.csv")

sim_res_1week_slow_late <- sim_SEIRAHD_Simulx_ME_2params(path_to_ind_params = paste0(getwd(), "/ind_params_lowb1.txt"), 
                                                    regressor_df_path = "ld1_reg_1week_slow_late.csv")
sim_res_2week_slow_late <- sim_SEIRAHD_Simulx_ME_2params(path_to_ind_params = paste0(getwd(), "/ind_params_lowb1.txt"), 
                                                    regressor_df_path = "ld1_reg_2week_slow_late.csv")
  
sim_res_slow_NPI <- list(sim_res_1week_slow, sim_res_2week_slow, 
                         sim_res_1week_slow_late, sim_res_2week_slow_late)

sim_res_slow_names <- c("1week_slow", "2week_slow", "1week_slow_late", "2week_slow_late")
  
for(j in 1:4){
  monolix_SEIRAHD <- sim_res_slow_NPI[[j]] %>%
    left_join(list_reg_files_slow_NPI[[j]] %>% select(id, time, lockdown1, BG1), 
              by = c("id", "time")) %>%
    rename(dept_id = id, day = time) %>%
    left_join(popsize_df, by = "dept_id") %>%
    group_by(dept_id) %>%
    mutate(initH = PrevH_ME[day == 1]) %>%
    ungroup() %>%
    select(dept_id, day, IncI_ME, IncH_ME, PrevH_ME, IncD_ME, initH, lockdown1, BG1) %>%
    pivot_longer(c(IncI_ME, IncH_ME, PrevH_ME, IncD_ME), 
                 names_to = "obs_id", values_to = "obs") %>%
    mutate(obs_id = case_when(obs_id == "IncH_ME" ~ 1, 
                              obs_id == "PrevH_ME" ~ 2,
                              obs_id == "IncI_ME" ~ 3, 
                              obs_id == "IncD_ME" ~ 4)) %>%
    relocate(obs, .after = day) %>%
    relocate(obs_id, .after = obs) %>%
    mutate(obs = ifelse(obs < 0, 0, obs))
  
  
  write.table(monolix_SEIRAHD, file = paste0("data_sim_SEIRAHD_2params_", sim_res_slow_names[j], ".txt"), 
              row.names = FALSE, sep = ",")
  
  monolix_SEIR <- monolix_SEIRAHD %>%
    filter(obs_id == 3) %>%
    rename(IncI = obs) %>%
    select(-c(obs_id, initH))
  
  write.table(monolix_SEIR, file = paste0("data_sim_SEIR_Simulx_2params_", sim_res_slow_names[j], ".txt"), 
              sep = ",", row.names = FALSE)
}

save(sim_res_slow_NPI, file = "sim_res_slow_NPI.RData")


for(j in 1:4){
  data <- read.table(paste0("data_sim_SEIR_Simulx_2params_", sim_res_slow_names[j], ".txt"), 
                     sep = ",", header = TRUE)
  
  plot <- ggplot(data, aes(x = day, y = IncI, group = dept_id)) + 
    geom_line() +
    labs(title = paste("Sim ", sim_res_slow_names[j]))
  
  
  print(plot)
}
