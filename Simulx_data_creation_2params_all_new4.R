library(tidyverse)
library(colorspace)
library(lixoftConnectors)
initializeLixoftConnectors(software = "simulx")

setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params")
dir <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims"


#### Simulx data creation ####
source(paste0(dir, "/SEIRAHD_Simulx_function_2params.R"))
source(paste0(dir, "/useful_functions.R"))

# load initial parameters
load("init_list.RData")


# change files with initial parameters
for(j in 1:100){
  indprms <- read.table(paste0(getwd(), "/ind_params/ind_2params_", j, ".txt"), 
                        header = TRUE, sep = " ") %>%
    mutate(initE_quantiles = cut(log_initE, breaks = 4, labels = 1:4), 
           log_initE = case_when(initE_quantiles == 1 ~ log_initE, 
                                 initE_quantiles == 2 ~ log_initE/2, 
                                 initE_quantiles == 3 ~ log_initE/3, 
                                 initE_quantiles == 4 ~ log_initE/4)) %>%
    select(-c(beta_w, beta_school, initE_quantiles))
  
  write.table(indprms, file = paste0(getwd(), "/ind_params/ind_2params_new4_", j, ".txt"), 
              row.names = F, sep = " ")
}

# change regressors (initH only)
later_reg_part <- data.frame(id = rep(1:94, each = 29), 
                             time = rep(122:150, 94), 
                             lockdown1 = 0, 
                             BG1 = 1)


reg_file <- read.csv("ld1_reg_df_2params.csv") %>%
  mutate(initH = initH/5) %>%
  group_by(id) %>%
  mutate(lockdown1 = lag(lockdown1, 29, default = 0), 
         BG1 = lag(BG1, 29, default = 0)) %>%
  ungroup() %>%
  bind_rows(., later_reg_part) %>%
  arrange(id, time) %>%
  tidyr::fill(initH, popsize)
  
write.csv(reg_file, file = "ld1_reg_df_2params_new4.csv", row.names = FALSE)


#### Simulx simulations new 4 ####
project.file <- "sim_SEIRAHD_2params_init_est.mlxtran"
importMonolixProject(project.file)

sim_res_Simulx_2params_new4_list <- list()

for(j in 1:100){
  path_to_ind_params <- paste0(getwd(), "/ind_params/ind_2params_new4_", j, ".txt")
  
  sim_res <- sim_SEIRAHD_Simulx_ME_2params(path_to_ind_params = path_to_ind_params, 
                                           regressor_df_path = "ld1_reg_df_2params_new4.csv",
                                           tmax = 150)
  
  sim_res_Simulx_2params_new4_list[[j]] <- sim_res
  
  monolix_SEIRAHD <- monolix_data_creation_ME(simulation_results = sim_res,
                                              popsize_df = popsize_df, 
                                              start_date = as.Date("2020-02-01"), 
                                              ld1_start = 45, ld1_end = 99, BG1_start = 100)
  
  write.table(monolix_SEIRAHD, file = paste0("data_sim_SEIRAHD_Simulx_2params_new4_ME", j, ".txt"), 
              row.names = FALSE, sep = ",")
  
  monolix_SEIR <- monolix_SEIRAHD %>%
    filter(obs_id == 3) %>%
    rename(IncI = obs) %>%
    select(-c(obs_id, initH))
  
  write.table(monolix_SEIR, file = paste0("data_sim_SEIR_Simulx_2params_new4_ME", j, ".txt"), 
              sep = ",", row.names = FALSE)
}

save(sim_res_Simulx_2params_new4_list, file = "sim_res_Simulx_2params_new4_list.RData")

### visualization ####
for(j in 1:100){
  data <- read.table(paste0("data_sim_SEIRAHD_Simulx_2params_new4_ME", j, ".txt"), 
                     sep = ",", header = TRUE)
  
  plot <- ggplot(data %>% filter(obs_id == 3), aes(x = day, y = obs*popsize/10^4, group = dept_id)) + 
    geom_line() 
  +
    scale_y_continuous(limits = c(0, 10))
  
  print(plot)
}


#### shorten data sets ####
for(j in 1:100){
  data <- read.table(paste0("data_sim_SEIRAHD_Simulx_2params_new4_ME", j, ".txt"), 
                     sep = ",", header = TRUE)
  
  data_short <- data %>%
    filter(day > 29) %>%
    mutate(day = day -29)
  
  write.table(data_short, file = paste0("data_sim_SEIRAHD_Simulx_2params_new4_short_ME", j, ".txt"), 
              sep = ",", row.names = FALSE)
  
  
  data_SEIR_short <- data_short %>%
    filter(obs_id == 3) %>%
    rename(IncI = obs) %>%
    select(-c(obs_id, initH))
  
  write.table(data_SEIR_short, file = paste0("data_sim_SEIR_Simulx_2params_new4_short_ME", j, ".txt"), 
              sep = ",", row.names = FALSE)
  
}
