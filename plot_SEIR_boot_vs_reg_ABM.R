library(tidyverse)
library(colorspace)
library(magrittr)
library(kableExtra)


# source("~/PhD/COVID_France/Dropbox_iris_covid/departement/Données_SPF/Data/data_functions.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/useful_functions.R")


setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots")

dir <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/sim_2params_regs"
dir3 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/boot_sim_2params_ABM5"
dir4 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/boot_sim_2params_ABM6"



#### ABM 5 ####
# load data
load(paste0(dir3, "/ABM5_hybrid_list.RData"))
load(paste0(dir3, "/ABM5_rm_list.RData"))

load(paste0(dir, "/reg_res_I_2params_all_ABM_rm5_df.RData"))
load(paste0(dir, "/reg_res_I_2params_all_ABM_hybrid5_df.RData"))



# bring data into right format and assemble
boot_df_ABM <- bootstrap_summary(ABM5_rm_list) %>% 
  mutate(model = "SEIR ABM rm") %>%
  bind_rows(bootstrap_summary(ABM5_hybrid_list) %>%
              mutate(model = "SEIR ABM hybrid")) %>%
  rename(mean_est = mean_est2, CI_LL = CI_LL2, CI_UL = CI_UL2)


reg_res_2params_ABM <- reg_res_I_2params_all_ABM_hybrid5_df %>%
  mutate(model = "Regression model ABM hybrid") %>%
  bind_rows(reg_res_I_2params_all_ABM_rm5_df %>% mutate(model = "Regression model ABM random mixing")) %>%
  reg_summary()

comp_df_ABM <- boot_df_ABM %>%
  bind_rows(reg_res_2params_ABM %>% rename(sim_rep = rep, mean_est = value))

# plot
pg_palette <- diverging_hcl("Purple-Green", n = 20)
comp_cols <- c(pg_palette[c(18, 14)], "black", pg_palette[c(3, 6)])

ggplot(comp_df_ABM, aes(x = sim_rep, y = mean_est, ymin = CI_LL, ymax = CI_UL, col = model)) + 
  geom_pointrange(position = position_dodge(width = 0.8)) +
  geom_line(aes(y = true_value), col = "black", linetype = "dashed", linewidth = 0.8) + 
  scale_x_continuous(breaks = c(1, seq(10, 100, 10)), expand = c(0.01, 0.01)) + 
  facet_wrap(~parameter, scales = "free_y", nrow = 2) +
  labs(title = "Datasets simulated with ABM", 
       x = "Simulation data set", y = "Coefficient value", col = "") + 
  scale_color_manual(values = comp_cols[-3]) +
  theme_bw() +
  theme(plot.title = element_text(family = "serif", size = 20), 
        axis.title = element_text(family = "serif", size = 15), 
        axis.text.x = element_text(family = "serif", size = 13), 
        axis.text.y = element_text(family = "serif", size = 13), 
        legend.text = element_text(family = "serif", size = 14.5))

ggsave("SEIR_vs_reg_ABM.jpeg", dpi = 400, width = 14, height = 8.5)

# bring into right format for metrics table
metric_df_ABM <- comp_df_ABM %>%
  select(parameter, model, perc_CI_covers, mean_bias, mean_rel_bias) %>%
  unique() %>%
  rename(coverage = perc_CI_covers, bias = mean_bias, relbias = mean_rel_bias) %>%
  pivot_longer(cols = c(coverage, bias, relbias), 
               names_to = "metric", 
               values_to = "value") %>%
  pivot_wider(names_from = model, values_from = value) %>%
  mutate(metric = factor(metric, levels = c("coverage", "bias", "relbias"))) %>%
  arrange(parameter, metric) %>%
  mutate(metric = case_when(metric == "bias" ~ "absolute bias", 
                            metric == "coverage" ~ "CI coverage (%)", 
                            metric == "relbias" ~ "relative bias (%)"))

# print metrics table
metric_df_ABM  %>% 
  ungroup() %>%
  select(-parameter) %>%
  kable(digits = 2, format = "html", table.attr = "style='width:45%;'") %>%
  kable_styling(bootstrap_options = "striped") %>%
  pack_rows("NPI 1", 1, 3) %>%
  pack_rows("NPI 2", 4, 6)




#### ABM 6 ####
# load data
load(paste0(dir4, "/ABM6_hybrid_list.RData"))
load(paste0(dir4, "/ABM6_rm_list.RData"))

load(paste0(dir, "/reg_res_I_2params_all_ABM_rm6_df.RData"))
load(paste0(dir, "/reg_res_I_2params_all_ABM_hybrid6_df.RData"))



# bring data into right format and assemble
boot_df_ABM6 <- bootstrap_summary(ABM6_rm_list, true_val_NPI2 = -0.8) %>% 
  mutate(model = "SEIR ABM rm") %>%
  bind_rows(bootstrap_summary(ABM6_hybrid_list, true_val_NPI2 = -0.8) %>%
              mutate(model = "SEIR ABM hybrid")) %>%
  rename(mean_est = mean_est2, CI_LL = CI_LL2, CI_UL = CI_UL2)


reg_res_2params_ABM6 <- reg_res_I_2params_all_ABM_hybrid6_df %>%
  mutate(model = "Regression model ABM hybrid") %>%
  bind_rows(reg_res_I_2params_all_ABM_rm6_df %>% mutate(model = "Regression model ABM random mixing")) %>%
  reg_summary(true_val_NPI2 = -0.8)

comp_df_ABM6 <- boot_df_ABM6 %>%
  bind_rows(reg_res_2params_ABM6 %>% rename(sim_rep = rep, mean_est = value))

# plot
pg_palette <- diverging_hcl("Purple-Green", n = 20)
comp_cols <- c(pg_palette[c(18, 14)], "black", pg_palette[c(3, 6)])

ggplot(comp_df_ABM6, aes(x = sim_rep, y = mean_est, ymin = CI_LL, ymax = CI_UL, col = model)) + 
  geom_pointrange(position = position_dodge(width = 0.8)) +
  geom_line(aes(y = true_value), col = "black", linetype = "dashed", linewidth = 0.8) + 
  scale_x_continuous(breaks = c(1, seq(10, 100, 10)), expand = c(0.01, 0.01)) + 
  facet_wrap(~parameter, scales = "free_y", nrow = 2) +
  labs(title = "Datasets simulated with ABM", 
       x = "Simulation data set", y = "Coefficient value", col = "") + 
  scale_color_manual(values = comp_cols[-3]) +
  theme_bw() +
  theme(plot.title = element_text(family = "serif", size = 20), 
        axis.title = element_text(family = "serif", size = 15), 
        axis.text.x = element_text(family = "serif", size = 13), 
        axis.text.y = element_text(family = "serif", size = 13), 
        legend.text = element_text(family = "serif", size = 14.5))

ggsave("SEIR_vs_reg_ABM6.jpeg", dpi = 400, width = 14, height = 8.5)

# bring into right format for metrics table
metric_df_ABM6 <- comp_df_ABM6 %>%
  select(parameter, model, perc_CI_covers, mean_bias, mean_rel_bias) %>%
  unique() %>%
  rename(coverage = perc_CI_covers, bias = mean_bias, relbias = mean_rel_bias) %>%
  pivot_longer(cols = c(coverage, bias, relbias), 
               names_to = "metric", 
               values_to = "value") %>%
  pivot_wider(names_from = model, values_from = value) %>%
  mutate(metric = factor(metric, levels = c("coverage", "bias", "relbias"))) %>%
  arrange(parameter, metric) %>%
  mutate(metric = case_when(metric == "bias" ~ "absolute bias", 
                            metric == "coverage" ~ "CI coverage (%)", 
                            metric == "relbias" ~ "relative bias (%)"))

# print metrics table
metric_df_ABM6  %>% 
  ungroup() %>%
  select(-parameter) %>%
  kable(digits = 2, format = "html", table.attr = "style='width:45%;'") %>%
  kable_styling(bootstrap_options = "striped") %>%
  pack_rows("NPI 1", 1, 3) %>%
  pack_rows("NPI 2", 4, 6)


