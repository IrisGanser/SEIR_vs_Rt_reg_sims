library(tidyverse)
library(colorspace)
library(lixoftConnectors)
initializeLixoftConnectors(software = "simulx")

setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params")


#### Simulx data creation ####
source("SEIRAHD_Simulx_function_2params.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/useful_functions.R")

# load initial parameters
load("init_list.RData")


# load and write regressors into file
popsize_df <- read.csv("popsize_df.csv") %>%
  mutate(dept_id = ifelse(dept_id < 20, dept_id, dept_id - 1))

ld1_reg_df_2params <- data.frame(id = rep(1:94, each = 121), 
                                 time = rep(1:121, 94), 
                                 initH = rep(init_list[[1]]$H, each = 121)) %>%
  mutate(lockdown1 = ifelse(between(time, 16, 70), 1, 0),
         BG1 = ifelse(time > 70, 1, 0), 
         popsize = 10000)

write.csv(ld1_reg_df_2params, "ld1_reg_df_2params.csv", row.names = F)


project.file <- "sim_SEIRAHD_2params_init_est.mlxtran"
importMonolixProject(project.file)

sim_res_Simulx_2params_new2_list <- list()

for(j in 1:100){
  path_to_ind_params <- paste0(getwd(), "/ind_params/ind_2params_", j, ".txt")
  
  sim_res <- sim_SEIRAHD_Simulx_ME_2params(path_to_ind_params = path_to_ind_params, 
                                           regressor_df_path = "ld1_reg_df_2params.csv")
  
  sim_res_Simulx_2params_new2_list[[j]] <- sim_res
  
  monolix_SEIRAHD <- monolix_data_creation_ME(simulation_results = sim_res,
                                              popsize_df = popsize_df)
  
  write.table(monolix_SEIRAHD, file = paste0("data_sim_SEIRAHD_Simulx_2params_new2_ME", j, ".txt"), 
              row.names = FALSE, sep = ",")
  
  monolix_SEIR <- monolix_SEIRAHD %>%
    filter(obs_id == 3) %>%
    rename(IncI = obs) %>%
    select(-c(obs_id, initH))
  
  write.table(monolix_SEIR, file = paste0("data_sim_SEIR_Simulx_2params_new2_ME", j, ".txt"), 
              sep = ",", row.names = FALSE)
}

save(sim_res_Simulx_2params_new2_list, file = "sim_res_Simulx_2params_new2_list.RData")

### visualization
for(j in 1:100){
  data <- read.table(paste0("data_sim_SEIRAHD_Simulx_2params_new2_ME", j, ".txt"), 
                     sep = ",", header = TRUE)
  
  plot <- ggplot(data %>% filter(obs_id == 3), aes(x = day, y = obs, group = dept_id)) + 
    geom_line()
  
  print(plot)
}
