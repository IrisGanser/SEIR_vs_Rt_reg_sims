library(tidyverse)
library(colorspace)
library(magrittr)
library(kableExtra)
library(patchwork)


source("~/PhD/COVID_France/Dropbox_iris_covid/departement/Données_SPF/Data/data_functions.R")
source("~/PhD/COVID_France/SEIR_vs_Rt_sims/useful_functions.R")


setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots")


dir5 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once"
dir6 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once2"
dir8 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once4"
dir7 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"


# load data
data_SEIRAHD <- read.table(paste0(dir7, "/data_sim_SEIRAHD_Simulx_2params_new2_ME1.txt"), 
                           header = TRUE, sep = ",")
data_ABM_rm_cov <- read.csv(paste0(dir8, "/data_covasim_rm5_Rt_1.csv")) %>%
  filter(day > 15) %>%
  mutate(day = day - 15)
data_ABM_hybrid_cov <- read.csv(paste0(dir8, "/data_covasim_hybrid5_Rt_1.csv")) %>%
  filter(day > 15) %>%
  mutate(day = day - 15)


# plot time series of cases and hospitalization and deaths (SEIRAHD only)
rect_cols <- sequential_hcl(5, palette = "BluYl")

p1 <- ggplot(data_SEIRAHD %>% filter(obs_id == 3), 
       aes(x = day, y = obs*popsize/10^4, group = dept_id)) +  
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = c(43, 96), y = Inf, label = c("NPI 1", "NPI 2"), 
           hjust = 0.5, vjust = 1, size = 4.5, fontface = 2, family = "serif") + 
  geom_line() + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 120, 10)) + 
  labs(title = "Simulx Cases", 
       x = "Day", y = "Cases") +
  theme_bw()  +
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12))

p1

p2 <- ggplot(data_SEIRAHD %>% filter(obs_id == 1), 
       aes(x = day, y = obs*popsize/10^4, group = dept_id)) +  
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = c(43, 96), y = Inf, label = c("NPI 1", "NPI 2"), 
           hjust = 0.5, vjust = 1, size = 4.5, fontface = 2, family = "serif") + 
  geom_line() + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 120, 10)) + 
  labs(title = "Simulx Hospitalizations", 
       x = "Day", y = "Hospital admissions") +
  theme_bw()  +
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12))

p2


p3 <- ggplot(data_ABM_rm_cov, aes(x = day, y = IncI*popsize/10^4, group = dept_id)) +  
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = c(43, 96), y = Inf, label = c("NPI 1", "NPI 2"), 
           hjust = 0.5, vjust = 1, size = 4.5, fontface = 2, family = "serif") + 
  geom_line() + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 120, 10)) + 
  scale_y_continuous(limits = c(0, 5000)) +
  labs(title = "Random mixing ABM", 
       x = "Day", y = "Cases") +
  theme_bw() +
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12))

p3

p4 <- ggplot(data_ABM_hybrid_cov, aes(x = day, y = IncI*popsize/10^4, group = dept_id)) +  
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = c(43, 96), y = Inf, label = c("NPI 1", "NPI 2"), 
           hjust = 0.5, vjust = 1, size = 4.5, fontface = 2, family = "serif") + 
  geom_line()  + 
  scale_y_continuous(limits = c(0, 5000)) + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 120, 10)) + 
  labs(title = "Multi-layer ABM", 
       x = "Day", y = "Cases") +
  theme_bw() +
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12))

p4

(p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 18, family = "serif", face = "bold", hjust = 0, vjust = 0))

ggsave("Data generation.jpeg", dpi = 400, width = 16, height = 9.5)

