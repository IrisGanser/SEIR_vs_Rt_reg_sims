library(tidyverse)
library(colorspace)
library(magrittr)
library(kableExtra)
library(patchwork)


source("~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIR_vs_Rt_reg_sims/useful_functions.R")


setwd("~/PhD/COVID_France/SEIR_vs_Rt_sims/plots")

dir1 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/SEIRAHD_Simulx_data_creation_2params"
dir5 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once4"
dir6 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once6"
dir7 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once7"
dir10 <- "~/PhD/COVID_France/SEIR_vs_Rt_sims/ABM_2params_all_at_once10"


# load data
data_SEIRAHD2 <- read.table(paste0(dir1, "/data_sim_SEIRAHD_Simulx_2params_new2_ME1.txt"), 
                           header = TRUE, sep = ",")
data_SEIRAHD3 <- read.table(paste0(dir1, "/data_sim_SEIRAHD_Simulx_2params_new3_ME1.txt"), 
                            header = TRUE, sep = ",")
data_SEIRAHD4 <- read.table(paste0(dir1, "/data_sim_SEIRAHD_Simulx_2params_new4_ME1.txt"), 
                            header = TRUE, sep = ",")

data_ABM_rm_cov5 <- read.csv(paste0(dir5, "/data_covasim_rm5_Rt_1.csv")) %>%
  filter(day > 15) %>%
  mutate(day = day - 15)
data_ABM_hybrid_cov5 <- read.csv(paste0(dir5, "/data_covasim_hybrid5_Rt_1.csv")) %>%
  filter(day > 15) %>%
  mutate(day = day - 15)

data_ABM_rm_cov6 <- read.csv(paste0(dir6, "/data_covasim_rm6_Rt_4.csv")) %>%
  filter(day > 15) %>%
  mutate(day = day - 15)
data_ABM_hybrid_cov6 <- read.csv(paste0(dir6, "/data_covasim_hybrid6_Rt_1.csv")) %>%
  filter(day > 15) %>%
  mutate(day = day - 15)

data_ABM_rm_cov7 <- read.csv(paste0(dir7, "/data_covasim_rm7_Rt_1.csv")) %>%
  filter(day > 29) %>%
  mutate(day = day - 29)
data_ABM_hybrid_cov7 <- read.csv(paste0(dir7, "/data_covasim_hybrid7_Rt_1.csv")) %>%
  filter(day > 29) %>%
  mutate(day = day - 29)

data_ABM_rm_cov10 <- read.csv(paste0(dir10, "/data_covasim_rm10_Rt_1.csv")) %>%
  filter(day > 29) %>%
  mutate(day = day - 29)
data_ABM_hybrid_cov10 <- read.csv(paste0(dir10, "/data_covasim_hybrid10_Rt_1.csv")) %>%
  filter(day > 29) %>%
  mutate(day = day - 29)

# plot time series of cases and hospitalization and deaths (SEIRAHD only)
rect_cols <- sequential_hcl(5, palette = "BluYl")

p1 <- ggplot(data_SEIRAHD2 %>% filter(obs_id == 3), 
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

p2 <- ggplot(data_SEIRAHD2 %>% filter(obs_id == 1), 
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


p3 <- ggplot(data_ABM_rm_cov5, aes(x = day, y = IncI*popsize/10^4, group = dept_id)) +  
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

p4 <- ggplot(data_ABM_hybrid_cov5, aes(x = day, y = IncI*popsize/10^4, group = dept_id)) +  
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



# plot with higher NPI 2 data

p5 <- ggplot(data_SEIRAHD3 %>% filter(obs_id == 3), 
             aes(x = day, y = obs*popsize/10^4, group = dept_id)) +  
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = c(43, 96), y = Inf, label = c("NPI 1", "NPI 2"), 
           hjust = 0.5, vjust = 1, size = 4.5, fontface = 2, family = "serif") + 
  geom_line() + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 120, 10)) + 
  labs(title = "SEIRAHD model cases", 
       x = "Day", y = "Cases") +
  theme_bw()  +
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12))

p5

p6 <- ggplot(data_SEIRAHD3 %>% filter(obs_id == 1), 
             aes(x = day, y = obs*popsize/10^4, group = dept_id)) +  
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = c(43, 96), y = Inf, label = c("NPI 1", "NPI 2"), 
           hjust = 0.5, vjust = 1, size = 4.5, fontface = 2, family = "serif") + 
  geom_line() + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 120, 10)) + 
  labs(title = "SEIRAHD model hospitalizations", 
       x = "Day", y = "Hospital admissions") +
  theme_bw()  +
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12))

p6


p7 <- ggplot(data_ABM_rm_cov6, aes(x = day, y = IncI, group = dept_id)) +  
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = c(43, 96), y = Inf, label = c("NPI 1", "NPI 2"), 
           hjust = 0.5, vjust = 1, size = 4.5, fontface = 2, family = "serif") + 
  geom_line() + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 120, 10)) + 
  labs(title = "Random mixing ABM", 
       x = "Day", y = "Cases") +
  theme_bw() +
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12))

p7

p8 <- ggplot(data_ABM_hybrid_cov6, aes(x = day, y = IncI, group = dept_id)) +  
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = c(43, 96), y = Inf, label = c("NPI 1", "NPI 2"), 
           hjust = 0.5, vjust = 1, size = 4.5, fontface = 2, family = "serif") + 
  geom_line()  + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 120, 10)) + 
  labs(title = "Multi-layer ABM", 
       x = "Day", y = "Cases") +
  theme_bw() +
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12))

p8

(p5 + p6) / (p7 + p8) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 18, family = "serif", face = "bold", hjust = 0, vjust = 0))

ggsave("Data generation high NPI 2.jpeg", dpi = 400, width = 16, height = 9.5)



# longer pre-analysis period (ABM7)
p9 <- ggplot(data_ABM_rm_cov7, aes(x = day, y = IncI, group = dept_id)) +  
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = c(43, 96), y = Inf, label = c("NPI 1", "NPI 2"), 
           hjust = 0.5, vjust = 1, size = 4.5, fontface = 2, family = "serif") + 
  geom_line() + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 120, 10)) + 
  labs(title = "Random mixing ABM", 
       x = "Day", y = "Cases") +
  theme_bw() +
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12))

p9

p10 <- ggplot(data_ABM_hybrid_cov7, aes(x = day, y = IncI, group = dept_id)) +  
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = c(43, 96), y = Inf, label = c("NPI 1", "NPI 2"), 
           hjust = 0.5, vjust = 1, size = 4.5, fontface = 2, family = "serif") + 
  geom_line()  + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 120, 10)) + 
  labs(title = "Multi-layer ABM", 
       x = "Day", y = "Cases") +
  theme_bw() +
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12))

p10

p9 + p10
ggsave("Data generation ABM6.jpeg", dpi = 400, width = 12, height = 5)



### SEIRAHD 4 ####
p11 <- ggplot(data_SEIRAHD4 %>% filter(obs_id == 3), 
             aes(x = day, y = obs*popsize/10^4, group = dept_id)) +  
  annotate("rect", xmin = 45, xmax = 100, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 100, xmax = 151, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = c(73, 126), y = Inf, label = c("NPI 1", "NPI 2"), 
           hjust = 0.5, vjust = 1, size = 4.5, fontface = 2, family = "serif") + 
  geom_line() + 
  geom_vline(xintercept = 10, linetype = "dashed") + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 150, 10)) + 
  labs(title = "Simulx Cases", 
       x = "Day", y = "Cases") +
  theme_bw()  +
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12))

p11

p12 <- ggplot(data_SEIRAHD4 %>% filter(obs_id == 1), 
             aes(x = day, y = obs*popsize/10^4, group = dept_id)) +  
  annotate("rect", xmin = 45, xmax = 100, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 100, xmax = 151, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = c(73, 126), y = Inf, label = c("NPI 1", "NPI 2"), 
           hjust = 0.5, vjust = 1, size = 4.5, fontface = 2, family = "serif") + 
  geom_line() + 
  geom_vline(xintercept = 10, linetype = "dashed") + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 120, 10)) + 
  labs(title = "Simulx Hospitalizations", 
       x = "Day", y = "Hospital admissions") +
  theme_bw()  +
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12))

p12

p11+p12



#### ABM 10 ####
p13 <- ggplot(data_ABM_rm_cov10, aes(x = day, y = IncI, group = dept_id)) +  
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = c(43, 96), y = Inf, label = c("NPI 1", "NPI 2"), 
           hjust = 0.5, vjust = 1, size = 4.5, fontface = 2, family = "serif") + 
  geom_line() + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 120, 10)) + 
  labs(title = "Random mixing ABM", 
       x = "Day", y = "Cases") +
  theme_bw() +
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12))

p13

p14 <- ggplot(data_ABM_hybrid_cov10, aes(x = day, y = IncI, group = dept_id)) +  
  annotate("rect", xmin = 16, xmax = 71, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[3]) +
  annotate("rect", xmin = 71, xmax = 121, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = rect_cols[4]) +
  annotate("label", x = c(43, 96), y = Inf, label = c("NPI 1", "NPI 2"), 
           hjust = 0.5, vjust = 1, size = 4.5, fontface = 2, family = "serif") + 
  geom_line()  + 
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 120, 10)) + 
  labs(title = "Multi-layer ABM", 
       x = "Day", y = "Cases") +
  theme_bw() +
  theme(plot.title = element_text(family = "serif", size = 16), 
        axis.title = element_text(family = "serif", size = 13), 
        axis.text.x = element_text(family = "serif", size = 12), 
        axis.text.y = element_text(family = "serif", size = 12))

p14


(p5 + p6) / (p13 + p14) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 18, family = "serif", face = "bold", hjust = 0, vjust = 0))

ggsave("Data generation BG0.8 models.jpeg", dpi = 400, width = 16, height = 9.5)
