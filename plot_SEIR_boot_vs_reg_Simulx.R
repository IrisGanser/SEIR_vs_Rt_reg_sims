library(tidyverse)
library(colorspace)
library(magrittr)
library(kableExtra)


source("~/PhD/COVID_France/Dropbox_iris_covid/departement/Données_SPF/Data/data_functions.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/useful_functions.R")


setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots")

dir <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/sim_2params_regs"
dir1 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/boot_sim_2params_all_new2_Simulx_SEIR"
dir2 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/boot_sim_2params_all_new2_Simulx_SEIRAHD"


# load data
load(paste0(dir1, "/SEIR_boot_list_2params_new2_Simulx.RData"))
load(paste0(dir2, "/SEIRAHD_boot_list_2params_new2_Simulx.RData"))
load(paste0(dir, "/reg_res_I_2params_all_Simulx_df.RData"))
load(paste0(dir, "/reg_res_H_2params_all_Simulx_df.RData"))


# bring data into right format and assemble
SEIR_boot_df_comp <- bootstrap_summary(SEIR_boot_list_all) 
SEIRAHD_boot_df_comp <- bootstrap_summary(SEIRAHD_boot_list_all) 

reg_res_2params_Simulx_new2 <- reg_res_I_2params_all_Simulx_df %>%
  mutate(model = "Regression model IncI") %>%
  bind_rows(reg_res_H_2params_all_Simulx_df %>% mutate(model = "Regression model IncH")) %>%
  reg_summary()


comp_df_Simulx <- SEIR_boot_df_comp %>%
  mutate(model = "SEIR model") %>%
  bind_rows(SEIRAHD_boot_df_comp %>% mutate(model = "SEIRAHD model")) %>%
  bind_rows(reg_res_2params_Simulx_new2 %>% rename(sim_rep = rep, mean_est = value)) %>%
  mutate(parameter = ifelse(parameter == "Lockdown", "NPI 1", "NPI 2"))


# plot
br_palette <- diverging_hcl("Blue-Red", n = 20)
comp_cols <- c(br_palette[c(1, 5)], "black", br_palette[c(16, 20)])

ggplot(comp_df_Simulx, aes(x = sim_rep, y = mean_est, ymin = CI_LL, ymax = CI_UL, col = model)) + 
  geom_pointrange(position = position_dodge(width = 0.8)) +
  geom_line(aes(y = true_value), col = "black", linetype = "dashed", linewidth = 0.8) + 
  facet_wrap(~parameter, scales = "free_y", nrow = 2) +
  labs(title = "95% CIs, SEIR-type models vs. 2-step regression models", 
       x = "Simulation data set", y = "Coefficient value", col = "") + 
  scale_color_manual(values = comp_cols[-3]) +
  theme_bw() +
  theme(plot.title = element_text(family = "serif", size = 20), 
        axis.title = element_text(family = "serif", size = 15), 
        axis.text.x = element_text(family = "serif", size = 13), 
        axis.text.y = element_text(family = "serif", size = 13), 
        legend.text = element_text(family = "serif", size = 14.5))

ggsave("SEIR_vs_reg_Simulx.jpeg", dpi = 400, width = 14, height = 8.5)


# metrics table
# bring into right format for metrics table
metric_df_Simulx <- comp_df_Simulx %>%
  select(parameter, model, perc_CI_covers, mean_bias, mean_rel_bias) %>%
  unique() %>%
  rename(coverage = perc_CI_covers, bias = mean_bias, relbias = mean_rel_bias) %>%
  pivot_longer(cols = c(coverage, bias, relbias), 
               names_to = "metric", 
               values_to = "value") %>%
  pivot_wider(names_from = model, values_from = value) %>%
  mutate(metric = factor(metric, levels = c("coverage", "bias", "relbias"))) %>%
  arrange (parameter, metric) %>%
  mutate(metric = case_when(metric == "bias" ~ "absolute bias", 
                            metric == "coverage" ~ "CI coverage (%)", 
                            metric == "relbias" ~ "relative bias (%)"))

# print metrics table
metric_df_Simulx %>% 
  ungroup() %>%
  select(-parameter) %>%
  kable(digits = 2, format = "html") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
  pack_rows("NPI 1", 1, 3) %>%
  pack_rows("NPI 2", 4, 6)

