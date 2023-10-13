library(tidyverse)
library(lme4)
library(EpiEstim)
library(colorspace)
library(parallel)
library(foreach)
library(doParallel)
library(magrittr)
library(RColorBrewer)

setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/Rt_trajectories")
source("~/PhD/COVID_France/Dropbox_iris_covid/departement/Données_SPF/Data/data_functions.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")


dir2 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"
dir3 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once"
dir4 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once2"


popsize_df <- read.csv(paste0(dir2, "/popsize_df.csv")) %>%
  mutate(dept_id = ifelse(dept_id < 20, dept_id, dept_id - 1))

load("ind_param_list.RData")


# EpiEstim bootstrap
data_list_Simulx <- loadRData(paste(dir2, "sim_res_Simulx_all_list2.RData", sep = "/"))

dataset1_Simulx <- data_list_Simulx[[1]] %>%
  left_join(popsize_df, by = c("id" = "dept_id")) %>%
  left_join(ind_param_list[[1]], by = "id") %>%
  mutate(Rt_SEIR = calc_Rt_SEIR(b1 = transmission, gamma = 1/5, S = S, N = 10000), 
         Rt_SEIRAHD = calc_Rt(b1 = transmission, S = S, Dq = 5, risk_hosp = 0.1, VE_I = 0, VE_H = 0), 
         IncI_unscaled = round(IncI*popsize/10^4), 
         IncH_unscaled = round(IncH*popsize/10^4), 
         lockdown1 = ifelse(between(time, 16, 70), 1, 0), 
         BG1 = ifelse(time > 70, 1, 0)) %>% 
  rename(dept_id = id, day = time)



cl <- makeCluster(10)
registerDoParallel(cl)

boot_EpiEstim <- list()

data_for_est = dataset1_Simulx
Inc_name = "IncI_unscaled"

seeds <- seq(101, 200)

for(j in 1:100){
  
  # resample individual departments from datasets
  set.seed(seeds[j])
  random_depts <- sample(1:94, 94, replace = TRUE)
  
  # Rt estimation
  Rt_list <- foreach(i = 1:94, .packages = c("tidyverse", "EpiEstim")) %dopar% {
    Inc_series <- data_for_est %>% 
      filter(dept_id == random_depts[i]) 
    
    Inc <- as.numeric(na.omit(Inc_series[[Inc_name]]*Inc_series$popsize/10^4))
    
    Rt_estim <- estimate_R(Inc,
                           method = "parametric_si",
                           config = make_config(list(
                             mean_si = 10.1,
                             std_si = 8.75)))$R
    Rt_estim <- Rt_estim %>%
      mutate(t = (t_end+t_start)/2,
             id = i, .before = t_start) %>%
      dplyr::select(t, id, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`) %>%
      rename(Rt = `Mean(R)`, CI_LL = `Quantile.0.025(R)`, CI_UL = `Quantile.0.975(R)`)
  }
  
  Rt_df <- do.call("rbind.data.frame", Rt_list) %>%
    rename(dept_id = id, day = t)
  
  # assemble dataset for regression
  reg_data <- data_for_est %>%
    dplyr::select(dept_id, day, lockdown1, BG1) %>%
    unique() %>%
    left_join(., Rt_df, by = c("dept_id", "day"))
    
    reg_res <- lmer(log(Rt) ~ lockdown1 + BG1 + (1|dept_id), data = reg_data)

  
  # get reg coefficients
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
    mutate(rep = j)
  
  boot_EpiEstim[[j]] <- coefs_Rt_reg
}

stopCluster(cl)

boot_EpiEstim_df <- do.call("rbind.data.frame", boot_EpiEstim)

boot_EpiEstim_res <- boot_EpiEstim_df %>%
  group_by(parameter) %>%
  summarize(CI_LL = quantile(value, probs = 0.025), 
            CI_UL = quantile(value, probs = 0.975), 
            value = mean(value))


### comparison: unbootstrapped run
cl <- makeCluster(8)
registerDoParallel(cl)

normal_reg_Simulx_I <- EpiEstim_reg_fun(data_for_est = dataset1_Simulx, 
                                     Inc_name = "IncI_unscaled", 
                                     rep_num = 1)
stopCluster(cl)


# compare both methods
boot_comp_df <- bind_rows(boot_EpiEstim_res %>% mutate(method = "bootstrap"), 
                          normal_reg_Simulx_I %>% mutate(method = "SE from reg only")) %>%
  mutate(true_value = factor(ifelse(parameter == "Lockdown 1", -1.45, -0.5)))

ggplot(boot_comp_df, aes(x = parameter, y = value, ymin = CI_LL, ymax = CI_UL, col = method)) + 
  geom_pointrange(position = position_dodge(width = 1)) + 
  #facet_wrap(~parameter, scales = "free_x") +
  geom_segment(aes(x = 0.5, xend = 1.5, y = -1.45, yend = -1.45), linetype = "dashed", col = "black") + 
  geom_segment(aes(x = 1.5, xend = 2.5, y = -0.5, yend = -0.5), linetype = "dashed", col = "black") + 
  labs(x = "", y = "coefficient value", 
       title = "Bootstrap vs. normal 2 step regresssion") + 
  theme_bw() + 
  scale_color_brewer(palette = "Dark2")
