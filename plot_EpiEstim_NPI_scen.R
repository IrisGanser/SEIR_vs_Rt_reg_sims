library(tidyverse)
library(colorspace)
library(magrittr)
library(kableExtra)
library(patchwork)
library(data.table)


source("~/PhD/COVID_France/Dropbox_iris_covid/departement/Données_SPF/Data/data_functions.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")


setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots")
dir9 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/Rt_trajectories"

load(paste0(dir9, "/Rt_NPI_scen_res_dfs.RData"))

res_reg_ld1_start_df %<>%
  mutate(parameter = ifelse(parameter == "Lockdown 1", "NPI 1", "NPI 2"))
res_reg_ld1_start_high_df %<>%
  mutate(parameter = ifelse(parameter == "Lockdown 1", "NPI 1", "NPI 2"))
res_reg_ld1_strength_df %<>%
  mutate(parameter = ifelse(parameter == "Lockdown 1", "NPI 1", "NPI 2"), 
         ld1_strength = 0-ld1_strength)

viridis_pal <- sequential_hcl(16, palette = "Viridis")

res_reg_ld1_start_df_copy <- data.table(res_reg_ld1_start_df) 
res_reg_ld1_start_high_df_copy <- data.table(res_reg_ld1_start_high_df) 

res_reg_ld1_start_df_copy[parameter == "NPI 1", y_min := -1.55]
res_reg_ld1_start_df_copy[parameter == "NPI 1", y_max := -0.8]
res_reg_ld1_start_df_copy[parameter == "NPI 2", y_min := -1.4]
res_reg_ld1_start_df_copy[parameter == "NPI 2", y_max := -0.49]

res_reg_ld1_start_high_df_copy[parameter == "NPI 1", y_min := -1.55]
res_reg_ld1_start_high_df_copy[parameter == "NPI 1", y_max := -0.8]
res_reg_ld1_start_high_df_copy[parameter == "NPI 2", y_min := -1.4]
res_reg_ld1_start_high_df_copy[parameter == "NPI 2", y_max := -0.49]



p1 <- ggplot(res_reg_ld1_start_df_copy, aes(ymin = CI_LL, ymax = CI_UL, x = ld1_start, 
                                 y = value, col = as.factor(ld1_start))) + 
  geom_pointrange(position = position_dodge(width = 1)) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~parameter, ncol = 1, scale = "free_y") +
  geom_line(aes(y = true_value), linetype = "dashed", col = "black", linewidth = 0.8) + 
  labs(title = "NPI start day (low transmission)", col = "NPI 1 start day",
       x = "NPI 1 start day", y = "Coefficient value") +
  theme_bw() + 
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max)) + 
  scale_color_brewer(palette = "Dark2") + 
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13),
        axis.text = element_text(family = "serif", size = 12), 
        legend.text = element_text(family = "serif", size = 12), 
        legend.title = element_text(family = "serif", size = 13), 
        strip.text = element_text(family = "serif", size = 13))


p2 <- ggplot(res_reg_ld1_start_high_df_copy, aes(ymin = CI_LL, ymax = CI_UL, x = ld1_start, 
                                      y = value, col = as.factor(ld1_start))) + 
  geom_pointrange(position = position_dodge(width = 1)) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~parameter, ncol = 1, scale = "free_y") +
  geom_line(aes(y = true_value), linetype = "dashed", col = "black", linewidth = 0.8) + 
  labs(title = "NPI start day (high transmission)", col = "NPI 1 start day",
       x = "NPI 1 start day", y = "Coefficient value") +
  theme_bw() +
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max)) + 
  scale_color_brewer(palette = "Dark2") + 
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13),
        axis.text = element_text(family = "serif", size = 12), 
        legend.text = element_text(family = "serif", size = 12), 
        legend.title = element_text(family = "serif", size = 13), 
        strip.text = element_text(family = "serif", size = 13))


p3 <- ggplot(res_reg_ld1_strength_df, aes(ymin = CI_LL, ymax = CI_UL, x = ld1_strength, 
                                    y = value, col = as.factor(ld1_strength))) + 
  geom_pointrange() + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(-2, -0.5, 0.1)) + 
  facet_wrap(~parameter, ncol = 1, scale = "free_y") +
  geom_segment(aes(y = true_value, yend = true_value, x = ld1_strength-0.05, xend = ld1_strength+0.05), 
               linetype = "dashed", col = "black", linewidth = 0.8) + 
  labs(title = "NPI strength", col = "NPI 1 coefficient",
       x = "NPI 1 strength", y = "Coefficient value") +
  theme_bw() +
  scale_color_manual(values = rev(viridis_pal), breaks = unique(res_reg_ld1_strength_df$ld1_strength)) + 
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13),
        axis.text = element_text(family = "serif", size = 12), 
        legend.text = element_text(family = "serif", size = 12), 
        legend.title = element_text(family = "serif", size = 13), 
        strip.text = element_text(family = "serif", size = 13))

plot1 <- p1 + p2 + plot_layout(guides = 'collect', tag_level = 'new')
plot1 / p3 + 
  plot_layout(heights = c(1, 1.3)) + 
  plot_annotation(tag_levels = c('A', '1')) & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 16, family = "serif", face = "bold", hjust = 0, vjust = 0))

ggsave("EpiEstim_NPI_scen.jpeg", dpi = 300, width = 14, height = 9)



p4 <- ggplot(res_EE_ld1_start_df %>% filter(dept_id %in% c(1, 4, 9, 13)),
       aes(x = day, y = Rt, col = as.factor(ld1_start))) + 
  geom_line(linewidth = 0.8, aes(linetype = "EpiEstim Rt")) +
  geom_ribbon(aes(ymax = CI_UL, ymin = CI_LL, fill = as.factor(ld1_start)), alpha = 0.5) +
  geom_line(aes(y = Rt_SEIRAHD, linetype = "Real Rt"), linewidth = 0.8) + 
  facet_wrap(~dept_id) + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 120, 30)) + 
  labs(title = "NPI start day (low transmission)", linetype = "", x = "Day",
       y = expression(R[t]), col = "NPI 1 start day", fill = "NPI 1 start day") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")  +
  scale_linetype(labels = c(expression(EpiEstim~R[t]), expression(Real~R[t]))) + 
  guides(color = guide_legend(order = 1),
         fill = guide_legend(order = 1),
         linetype = guide_legend(order = 2)) +
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text = element_text(family = "serif", size = 11), 
        legend.text = element_text(family = "serif", size = 12), 
        legend.title = element_text(family = "serif", size = 12.5), 
        legend.text.align = 0,
        strip.text = element_text(family = "serif", size = 12))


p5 <- ggplot(res_EE_ld1_start_high_df %>% filter(dept_id %in% c(1, 4, 9, 13)),
       aes(x = day, y = Rt, col = as.factor(ld1_start))) + 
  geom_line(linewidth = 0.8, aes(linetype = "EpiEstim Rt")) +
  geom_ribbon(aes(ymax = CI_UL, ymin = CI_LL, fill = as.factor(ld1_start)), alpha = 0.5) +
  geom_line(aes(y = Rt_SEIRAHD, linetype = "Real Rt"), linewidth = 0.8) + 
  facet_wrap(~dept_id) + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 120, 30)) + 
  labs(title = "NPI start day (high transmission)", linetype = "", x = "Day",
       y = expression(R[t]), col = "NPI 1 start day", fill = "NPI 1 start day") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")  +
  scale_linetype(labels = c(expression(EpiEstim~R[t]), expression(Real~R[t]))) +
  guides(color = guide_legend(order = 1),
         fill = guide_legend(order = 1),
         linetype = guide_legend(order = 2)) +
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text = element_text(family = "serif", size = 11), 
        legend.title = element_text(family = "serif", size = 12.5), 
        legend.text = element_text(family = "serif", size = 12),  
        legend.text.align = 0,
        strip.text = element_text(family = "serif", size = 12)) 


res_EE_ld1_strength_df %<>%
  mutate(strength_cat = cut(ld1_strength, 4, 
                            labels = c("Low", "Moderate low", "Moderate high", "High")), 
         ld1_strength = 0-ld1_strength)


p6 <- ggplot(res_EE_ld1_strength_df %>% filter(dept_id %in% c(1, 4, 9, 13)),
       aes(x = day, y = Rt, col = as.factor(ld1_strength))) + 
  geom_line(aes(y = Rt_SEIRAHD, linetype = "Real Rt"), linewidth = 0.8) + 
  geom_line(linewidth = 0.8, aes(linetype = "EpiEstim Rt")) +
  geom_ribbon(aes(ymax = CI_UL, ymin = CI_LL, fill = as.factor(ld1_strength)), alpha = 0.5) +
  facet_grid(rows = vars(strength_cat), cols = vars(dept_id)) + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 120, 30)) + 
  labs(title = "NPI strength", x = "Day",
       y = expression(R[t]), col = "NPI 1 coefficient", fill = "NPI 1 coefficient") +
  theme_bw() +
  scale_color_manual(values = rev(viridis_pal), breaks = unique(res_EE_ld1_strength_df$ld1_strength)) +
  scale_fill_manual(values = rev(viridis_pal), breaks = unique(res_EE_ld1_strength_df$ld1_strength)) + 
  guides(linetype = "none") + 
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text = element_text(family = "serif", size = 11), 
        legend.title = element_text(family = "serif", size = 12.5), 
        legend.text = element_text(family = "serif", size = 12),  
        legend.text.align = 0,
        strip.text = element_text(family = "serif", size = 12))


plot4a <- p4 + p5 + plot_layout(guides = 'collect', tag_level = 'new')
plot4a / p6 + 
  plot_layout(heights = c(1, 1.3)) + 
  plot_annotation(tag_levels = c('A', '1')) & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 16, family = "serif", face = "bold", hjust = 0, vjust = 0))

ggsave("Rt_traj_NPI_scen.jpeg", dpi = 300, width = 14, height = 9)

