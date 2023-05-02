library(here)
library(broom)
library(ggrepel)
library(lme4)
library(lmerTest)
library(sjPlot)
library(RLRsim)
library(MuMIn)
library(piecewiseSEM)
library(tidyverse)
library(patchwork)
library(ggthemes)
library(DiagrammeR)

theme_set(theme_base() + 
            theme(plot.background = element_blank(),
                  legend.title = element_text(size = 10),
                  legend.text = element_text(size = 8)))

pal <- colorspace::darken(RColorBrewer::brewer.pal(n = 10, name = "Set3"), amount = .2)

# 1. IMPORT DATA SETS ------------------------------------------------------------
#source("analysis_wrapper.R")
metacom_var <- read_csv(here("data/L4_metacommunity_variability_analysis_results_2023-04-19.csv"))
local_var <- read_csv(here("data/L4_local_variability_analysis_results_2023-04-19.csv"))
env_var <- read_csv(here("data/lter_centroid_satdata.csv"))
data_list <- read_csv(here("data/L3_DATA_list.csv")) %>% 
  filter(l3_filename != "L3-mcr-fish-castorani.csv")

metacom_var <- data_list %>% 
  rename(dataset_file_name = l3_filename) %>% 
  left_join(metacom_var, by = "dataset_file_name") %>% 
  rename(lter_site = `LTER site`)

organism_group_key <- (metacom_var %>% 
    select(dataset_id, organism_group)) %>% distinct()

local_var <- data_list %>% 
  rename(dataset_file_name = l3_filename) %>% 
  left_join(local_var, by = "dataset_file_name") %>% 
  rename(lter_site = `LTER site`)

unique(local_var$dataset_id) #33
unique(metacom_var$dataset_id) #33

# set up data frames

metacom_divstab_comp_dat <- metacom_var %>% select(dataset_id, lter_site, organism_group, variability_type, standardization_method, metric, metric_value) %>% 
  filter(standardization_method == "h") %>% 
  pivot_wider(names_from = "metric", values_from = "metric_value") %>% 
  left_join(metacom_var %>% 
              filter(variability_type == "divpart_time_series") %>% 
              select(dataset_id, metric, metric_value) %>% 
              pivot_wider(names_from = "metric", values_from = "metric_value"), 
            by = "dataset_id")

metacom_divstab_agg_dat <- metacom_var %>% select(dataset_id, lter_site, organism_group, variability_type, metric, metric_value) %>% 
  filter(variability_type == "agg") %>% 
  pivot_wider(names_from = "metric", values_from = "metric_value") %>% 
  left_join(metacom_var %>% 
              filter(variability_type == "divpart_time_series") %>% 
              select(dataset_id, metric, metric_value) %>% 
              pivot_wider(names_from = "metric", values_from = "metric_value"), 
            by = "dataset_id")

sum(metacom_divstab_agg_dat$dataset_id %in% metacom_divstab_comp_dat$dataset_id)

local_var <- local_var %>% filter(dataset_id %in% metacom_divstab_agg_dat$dataset_id)

local_divstab_regs <-  local_var %>% 
  select(dataset_id, lter_site, SITE_ID, organism_group, metric, metric_value) %>% 
  pivot_wider(names_from = "metric", values_from = "metric_value")

unique(local_divstab_regs$dataset_id)

# 2. LOCAL DIVERSITY STABILITY -----------------------------------------------
# This shows patterns within metacommunities about local community richness and variability
# Then compares across all sites, to show that sometimes these local relationships are stronger than others


### mixed effects models
local_dataset_for_mods <- local_var %>% select(dataset_id, lter_site, SITE_ID, organism_group, metric, metric_value) %>% 
  pivot_wider(names_from = "metric", values_from = "metric_value") 
local_dataset_for_mods <- local_dataset_for_mods %>% group_by(dataset_id) %>% 
  mutate(alpha_div_centered = site_mean_alpha_div - mean(site_mean_alpha_div),
         alpha_div_scaled = scale(site_mean_alpha_div))



local_div_stab_comp_alpha_mod <- lm(BD ~ alpha_div_scaled, data = local_dataset_for_mods)
local_div_stab_comp_alpha_fit <- glance(local_div_stab_comp_alpha_mod)
summary(local_div_stab_comp_alpha_mod)
plot(local_div_stab_comp_alpha_mod)
(p_val <- as.character(round(local_div_stab_comp_alpha_fit$p.value,2))) # 0.16
(r2 <- as.character(round(local_div_stab_comp_alpha_fit$r.squared,2))) # 0.


local_div_stab_agg_alpha_mod <- lm(CV ~ alpha_div_scaled, data = local_dataset_for_mods)
local_div_stab_agg_alpha_fit <- glance(local_div_stab_agg_alpha_mod)
summary(local_div_stab_agg_alpha_mod)
plot(local_div_stab_agg_alpha_mod)
(p_val <- as.character(round(local_div_stab_agg_alpha_fit$p.value,2))) # 0.
(r2 <- as.character(round(local_div_stab_agg_alpha_fit$r.squared,2))) # 0.09




#local_comp_mod_lmm <- lmer(BD ~ alpha_div_scaled + (alpha_div_scaled||dataset_id), data = local_dataset_for_mods)
local_comp_mod_lmm <- lmer(BD ~ alpha_div_scaled + (1 + alpha_div_scaled|dataset_id), data = local_dataset_for_mods)

summary(local_comp_mod_lmm)
plot(local_comp_mod_lmm)
coef(local_comp_mod_lmm)
confint(local_comp_mod_lmm, method = "boot")


# local_agg_mod_lmm <- lmer(CV ~ alpha_div_scaled + (alpha_div_scaled||dataset_id), data = local_dataset_for_mods)
# AIC(local_agg_mod_lmm)
local_agg_mod_lmm <- lmer(CV ~ alpha_div_scaled + (1 + alpha_div_scaled|dataset_id), data = local_dataset_for_mods)
AIC(local_agg_mod_lmm)
summary(local_agg_mod_lmm)
plot(local_agg_mod_lmm)
coef(local_agg_mod_lmm)
confint(local_agg_mod_lmm, method = "boot")

(r2m_comp <- r.squaredGLMM(local_comp_mod_lmm)[1])
(r2c_comp <- r.squaredGLMM(local_comp_mod_lmm)[2])

(r2m_agg <- r.squaredGLMM(local_agg_mod_lmm)[1])
(r2c_agg <- r.squaredGLMM(local_agg_mod_lmm)[2])

### make figs

R2m_expression <- expression(paste(" ", R[m]^2 , "= ", 0.013))
R2c_expression <- expression(paste(" ", R[c]^2 , "= ", 0.727))

local_divstab_comp_fig <- local_dataset_for_mods %>% 
                          
  ggplot(aes(x = alpha_div_scaled, y = BD)) + 
  geom_point(mapping = aes(group = dataset_id, color = organism_group), alpha = 0.3) + 
  geom_smooth(mapping = aes(group= dataset_id, color = organism_group), method = "lm",size=0.5, se = F, show.legend = FALSE) + 
  geom_smooth(method = "lm", se = F, color = "black", linewidth = 2, linetype = "dashed") +
  labs(x = expression(paste("Mean ", alpha, "-diversity (z-score)")), y = expression(paste("Comp. ", alpha, "-variability (", BD[alpha]^h,")")), color = "Organism group") + 
  scale_color_manual(values = pal) +
  #scale_y_log10() +
  annotate("text", x = 2.5, y = 0.8, label = R2m_expression) +
  annotate("text", x = 2.5, y = 0.7, label = R2c_expression)

local_divstab_comp_fig


R2m_expression <- expression(paste(" ", R[m]^2 , "= ", 0.022))
R2c_expression <- expression(paste(" ", R[c]^2 , "= ", 0.744))

local_divstab_agg_fig <- local_dataset_for_mods %>% 
  
  ggplot(aes(x = alpha_div_scaled, y = CV)) + 
  geom_point(mapping = aes(group = dataset_id, color = organism_group), alpha = 0.3) + 
  geom_smooth(mapping = aes(group = dataset_id, color = organism_group), method = "lm", size =0.5, se = F, show.legend = FALSE) + 
  geom_smooth(method = "lm", se = F, color = "black", size = 1.5) +
  labs(x = expression(paste("Mean ", alpha, "-diversity (z-score)")), y = expression(paste("Agg. ", alpha, "-variability (CV)")), color = "Organism group") + 
  scale_color_manual(values = pal) +
  #scale_y_log10() +
  annotate("text", x = 2.5, y = 1.6, label = R2m_expression) +
  annotate("text", x = 2.5, y = 1.4, label = R2c_expression)
local_divstab_agg_fig


local_divstab_fig <- local_divstab_agg_fig + local_divstab_comp_fig + 
  plot_layout(ncol = 1, guides = "collect") + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
ggsave(filename = here("figs/local_divstab_fig.png"), plot = local_divstab_fig, dpi = 600, width = 6, height = 6*3/4*2, bg = "white")



# 3. REGIONAL DIVERSITY STABILITY ------------------------------------------------
# What happens when we look at the regional scale?
# Do we see richness-variability relationships with composition and aggregate at across metacommunities?

div_stab_comp_gamma_mod <- (lm(gamma_var ~ gamma_div_mean, data = metacom_divstab_comp_dat))
div_stab_comp_beta_mod <- (lm(phi_var ~ beta_div_mean, data = metacom_divstab_comp_dat))
div_stab_comp_alpha_mod <- (lm(alpha_var ~ alpha_div_mean, data = metacom_divstab_comp_dat))
div_stab_comp_alpha_gamma_mod <- (lm(gamma_var ~ alpha_div_mean, data = metacom_divstab_comp_dat))
div_stab_comp_beta_gamma_mod <- (lm(gamma_var ~ beta_div_mean, data = metacom_divstab_comp_dat))

summary(div_stab_comp_alpha_mod)
summary(div_stab_comp_beta_mod)
summary(div_stab_comp_gamma_mod)
summary(div_stab_comp_alpha_gamma_mod)
summary(div_stab_comp_beta_gamma_mod)

div_stab_agg_gamma_mod <- (lm(gamma_var ~ gamma_div_mean, data = metacom_divstab_agg_dat))
div_stab_agg_beta_mod <- (lm(phi_var ~ beta_div_mean, data = metacom_divstab_agg_dat))
div_stab_agg_alpha_mod <- (lm(alpha_var ~ alpha_div_mean, data = metacom_divstab_agg_dat))
div_stab_agg_alpha_gamma_mod <- (lm(gamma_var ~ alpha_div_mean, data = metacom_divstab_agg_dat))
div_stab_agg_beta_gamma_mod <- (lm(gamma_var ~ beta_div_mean, data = metacom_divstab_agg_dat))

summary(div_stab_agg_alpha_mod)
summary(div_stab_agg_beta_mod)
summary(div_stab_agg_gamma_mod)
summary(div_stab_agg_alpha_gamma_mod)
summary(div_stab_agg_beta_gamma_mod)

# regional compositional variability

div_stab_comp_gamma_fit <- glance(div_stab_comp_gamma_mod)
(p_val <- as.character(round(div_stab_comp_gamma_fit$p.value,2))) # 0.05
(r2 <- as.character(round(div_stab_comp_gamma_fit$r.squared,2))) # 0.17

(div_stab_gamma_h <- metacom_divstab_comp_dat %>% 
    ggplot(aes(x = gamma_div_mean, y = gamma_var, label = lter_site)) +
    stat_smooth(method = "lm", se = T, size = 1, color = "black") +
    geom_point(size = 2, alpha = 0.7, mapping = aes(color = organism_group)) +
    #geom_label_repel(size = 2) +
    labs(x = expression(paste("Mean ", gamma, "-diversity")),
         y = expression(paste("Comp. ", gamma, "-variability (", BD^h [gamma], ")")),
         color = "Organism group")  +
    #scale_x_log10() +
    scale_color_manual(values = pal, drop = FALSE) +
    annotate("text", x = 10, y = 0.045, label = bquote(atop(paste(R^2, "= 0.04", ),
                                                            "p = 0.15")))
  #ggsave("ESA_2019/figs/variability_alpha-gamma.png", width = 6, height = 4, units = "in", dpi = 600)
)


div_stab_comp_alpha_fit <- glance(div_stab_comp_alpha_mod)
(p_val <- as.character(round(div_stab_comp_alpha_fit$p.value,2))) # 0.86
(r2 <- as.character(round(div_stab_comp_alpha_fit$r.squared,2))) # 0.
(div_stab_alpha_h <- metacom_divstab_comp_dat %>% 
    ggplot(aes(x = alpha_div_mean, y = alpha_var, label = lter_site)) +
    #stat_smooth(method = "lm", se = T, size = 1, color = "black") +
    geom_point(size = 2, alpha = 0.7, mapping = aes(color = organism_group)) +
    #geom_label_repel(size = 2) +
    labs(x = expression(paste("Mean ", alpha, "-diversity")),
         y = expression(paste("Comp. ", alpha, "-variability")),
         color = "Organism group") +
    scale_color_manual(values = pal, drop = FALSE) +
    #scale_x_log10() +
    annotate("text", x = 30, y = 0.09, label = "p = 0.87")
  #ggsave("ESA_2019/figs/variability_alpha-gamma.png", width = 6, height = 4, units = "in", dpi = 600)
)


div_stab_alpha_gamma_fit <- glance(div_stab_comp_alpha_gamma_mod)
div_stab_alpha_gamma_fit$p.value # 0.023
div_stab_alpha_gamma_fit$r.squared # 0.232
(div_stab_alpha_gamma_h <- metacom_divstab_comp_dat %>% 
    ggplot(aes(x = alpha_div_mean, y = gamma_var, label = lter_site)) +
    stat_smooth(method = "lm", se = T, size = 1, color = "black") +
    geom_point(size = 2, alpha = 0.7, mapping = aes(color = organism_group)) +
    #geom_label_repel(size = 2) +
    labs(x = expression(paste("Mean ", alpha, "-diversity")),
         y = expression(paste("Comp. ", gamma, "-variability")),
         color = "Organism group") +
    scale_color_manual(values = pal, drop = FALSE) +
    #scale_x_log10()  +
    annotate("text", x = 3, y = 0.045, label = bquote(atop(paste(R^2, "= 0.27"),
                                                           "p = 0.004")))  #ggsave("ESA_2019/figs/variability_alpha-gamma.png", width = 6, height = 4, units = "in", dpi = 600)
)

div_stab_beta_gamma_fit <- glance(div_stab_comp_beta_gamma_mod)
div_stab_beta_gamma_fit$p.value # 0.2
div_stab_beta_gamma_fit$r.squared
(div_stab_beta_h <- metacom_divstab_comp_dat %>% 
    ggplot(aes(x = beta_div_mean, y = gamma_var, label = lter_site)) +
    #stat_smooth(method = "lm", se = T, size = 1, color = "black") +
    geom_point(size = 2, alpha = 0.7, mapping = aes(color = organism_group)) +
    scale_color_manual(values = pal, drop = FALSE) +
    #geom_label_repel(size = 2) +
    labs(x = expression(paste("Mean ", beta, "-diversity")),
         y = expression(paste("Comp. ", gamma, "-variability")),
         color = "Organism group") +
    annotate("text", x = 4.75, y = 0.045, label = bquote("p = 0.34"))
  #ggsave("ESA_2019/figs/variability_alpha-gamma.png", width = 6, height = 4, units = "in", dpi = 600)
)

# regional aggregate variability
div_stab_gamma_agg_fit <- glance(div_stab_agg_gamma_mod)
div_stab_gamma_agg_fit$p.value # 0.202
div_stab_gamma_agg_fit$r.squared # 0.08
(div_stab_gamma_agg <- metacom_divstab_agg_dat %>% 
    ggplot(aes(x = gamma_div_mean, y = gamma_var, label = lter_site)) +
    #stat_smooth(method = "lm", se = T, size = 1, color = "black") +
    geom_point(size = 2, alpha = 0.7, mapping = aes(color = organism_group)) +
    scale_color_manual(values = pal, drop = FALSE) +
    #geom_label_repel(size = 2) +
    labs(x = expression(paste("Mean ", gamma, "-diversity")),
         y = expression(paste("Agg. ", gamma, "-variability (", CV^2, ")")),
         color = "Organism group") +
    annotate("text", x = 65, y = 0.025, label = bquote("p = 0.22"))
  #ggsave("ESA_2019/figs/variability_alpha-gamma.png", width = 6, height = 4, units = "in", dpi = 600)
)

div_stab_alpha_agg_fit <- glance(div_stab_agg_alpha_mod)
div_stab_alpha_agg_fit$p.value # 0.17
div_stab_alpha_agg_fit$r.squared # 0.0927
(div_stab_alpha_agg <- metacom_divstab_agg_dat %>% 
    ggplot(aes(x = alpha_div_mean, y = alpha_var, label = lter_site)) +
    #stat_smooth(method = "lm", se = T, size = 1, color = "black") +
    geom_point(size = 2, alpha = 0.7, mapping = aes(color = organism_group)) +
    scale_color_manual(values = pal, drop = FALSE) +
    #geom_label_repel(size = 2) +
    labs(x = expression(paste("Mean ", alpha, "-diversity")),
         y = expression(paste("Agg. ", alpha, "-variability")),
         color = "Organism group")  +
    # scale_x_log10() +
    annotate("text", x = 30, y = .09, label = bquote("p = 0.16"))
  #ggsave("ESA_2019/figs/variability_alpha-gamma.png", width = 6, height = 4, units = "in", dpi = 600)
)

div_stab_alpha_gamma_agg_fit <- glance(div_stab_agg_alpha_gamma_mod)
div_stab_alpha_gamma_agg_fit$p.value # 0.927
div_stab_alpha_gamma_agg_fit$r.squared # 0
(div_stab_alpha_gamma_agg <- metacom_divstab_agg_dat %>% 
    ggplot(aes(x = alpha_div_mean, y = gamma_var, label = lter_site)) +
    #stat_smooth(method = "lm", se = T, size = 1, color = "black") +
    geom_point(size = 2, alpha = 0.7, mapping = aes(color = organism_group)) +
    scale_color_manual(values = pal, drop = FALSE) +
    #geom_label_repel(size = 2) +
    labs(x = expression(paste("Mean ", alpha, "-diversity")),
         y = expression(paste("Agg. ", gamma, "-variability")),
         color = "Organism group")  +
    #scale_x_log10() +
    annotate("text", x = 30, y = 0.025, label = "p = 0.48")
  #ggsave("ESA_2019/figs/variability_alpha-gamma.png", width = 6, height = 4, units = "in", dpi = 600)
)

div_stab_beta_gamma_agg_fit <- glance(div_stab_agg_beta_gamma_mod)
div_stab_beta_gamma_agg_fit$p.value # 0.0095
div_stab_beta_gamma_agg_fit$r.squared # 0.29
(div_stab_beta_agg <- metacom_divstab_agg_dat %>% 
    ggplot(aes(x = beta_div_mean, y = gamma_var, label = lter_site)) +
    stat_smooth(method = "lm", se = T, size = 1, color = "black") +
    geom_point(size = 2, alpha = 0.7, mapping = aes(color = organism_group)) +
    scale_color_manual(values = pal, drop = FALSE) +
    #geom_label_repel(size = 2) +
    labs(x = expression(paste("Mean ", beta, "-diversity")),
         y = expression(paste("Agg. ", gamma, "-variability")),
         color = "Organism group")  + 
    annotate("text", x = 4.75, y = 0.02, label = bquote(atop(paste(R^2, "= 0.25"),
                                                               "p = 0.006")))
  #ggsave("ESA_2019/figs/variability_alpha-gamma.png", width = 6, height = 4, units = "in", dpi = 600)
)

# combine figs
div_stab_multipanel_fig <- div_stab_gamma_agg + div_stab_gamma_h +
  div_stab_beta_agg + div_stab_beta_h +
  div_stab_alpha_gamma_agg + div_stab_alpha_gamma_h +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect", nrow = 3)
ggsave(filename = "figs/diversity_variability.png", plot = div_stab_multipanel_fig, width = 2.5*4, height = 3*3, units = "in", dpi = 1000, bg = "white")


summary(div_stab_agg_gamma_mod)
summary(div_stab_comp_gamma_mod)

# 4. COMPARE BETA DIVERSITY WITH PHI, SYNCHRONY
summary(div_stab_comp_beta_mod)
div_stab_phi_fit <- glance(div_stab_comp_beta_mod)
div_stab_phi_fit$p.value # 2 e-4
div_stab_phi_fit$r.squared # 0.49
(div_stab_phi_h <- metacom_divstab_comp_dat %>% 
    ggplot(aes(x = beta_div_mean, y = phi_var, label = lter_site)) +
    stat_smooth(method = "lm", se = T, size = 1, color = "black") +
    geom_point(size = 2, alpha = 0.7, mapping = aes(color = organism_group)) +
    scale_y_continuous(limits = c(0,1)) +
    scale_color_manual(values = pal, drop = FALSE) +
    #geom_label_repel(size = 2) +
    labs(x = expression(paste("Mean ", beta, "-diversity")),
         y = expression(paste("Comp. Synchrony (", BD^h[phi], ")")),
         color = "Organism group")  +
    annotate("text", x = 4.75, y = 0.9, label = bquote(atop(paste(R^2, "= 0.38"),
                                                            "p = 0.0003")))
  #ggsave("ESA_2019/figs/variability_alpha-gamma.png", width = 6, height = 4, units = "in", dpi = 600)
)

summary(div_stab_agg_beta_mod)
div_stab_phi_agg_fit <- glance(div_stab_agg_beta_mod)
div_stab_phi_agg_fit$p.value # 0.006
div_stab_phi_agg_fit$r.squared # 0.32
(div_stab_phi_agg <- metacom_divstab_agg_dat %>% 
    ggplot(aes(x = beta_div_mean, y = phi_var, label = lter_site)) +
    stat_smooth(method = "lm", se = T, size = 1, color = "black") +
    geom_point(size = 2, alpha = 0.7, mapping = aes(color = organism_group)) +
    scale_y_continuous(limits = c(0,1))+
    scale_color_manual(values = pal, drop = FALSE) +
    #geom_label_repel(size = 2) +
    labs(x = expression(paste("Mean ", beta, "-diversity")),
         y = expression(paste("Agg. Synchrony (",phi,")")),
         color = "Organism group") +
    annotate("text", x = 4.75, y = 0.9, label = bquote(atop(paste(R^2, "= 0.35"),
                                                            "p = 0.006")))
  #ggsave("ESA_2019/figs/variability_alpha-gamma.png", width = 6, height = 4, units = "in", dpi = 600)
)


# combine figs
div_stab_phi_beta_fig <- (div_stab_phi_agg) + div_stab_phi_h + 
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", nrow = 2)
ggsave(filename = "figs/beta_phi.png", plot = div_stab_phi_beta_fig, width = 1.5*4, height = 2*3, units = "in", dpi = 1000, bg = "white")

# Compare alpha-alpha, beta-beta-, gamma-gamma
div_stab_multiscale_partition_fig <- 
  (local_divstab_agg_fig + theme(legend.position = "none")) + (local_divstab_comp_fig + theme(legend.position = "none")) +
  div_stab_phi_agg + div_stab_phi_h +
  #div_stab_gamma_agg + div_stab_gamma_h +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect", nrow = 2)
ggsave(filename = "figs/diversity_variability_alpha_beta.png", plot = div_stab_multiscale_partition_fig, width = 2.7*4, height = 3.2*3, units = "in", dpi = 1000, bg = "white")



# 4. COMPOSITIONAL VS AGGREGATE VARIABILITY -----------------------------------
# 

comp_agg_stab <- left_join(metacom_divstab_agg_dat, metacom_divstab_comp_dat, by = c("dataset_id", "lter_site", "organism_group"), suffix = c("_agg", "_comp"))


comp_agg_fig <- na.omit(comp_agg_stab) %>% 
  ggplot(aes(label = lter_site, 
             color = organism_group, group = paste(lter_site, organism_group))) + 
  #geom_point(mapping = aes(x = alpha_var_rate_comp, y = alpha_var_rate_agg), alpha = 0.5, size = 3, shape = 22) +
  #geom_point(mapping = aes(x = gamma_var_rate_comp, y = gamma_var_rate_agg), alpha = 0.5, size = 3, shape = 19) +
  geom_segment(aes(x = alpha_var_rate_comp, xend = gamma_var_rate_comp, 
                   y = alpha_var_rate_agg, yend = gamma_var_rate_agg),
               alpha = 0.7, arrow = arrow(length = unit(.2, "cm"), type = "closed", angle = 15)) +
  #geom_text_repel(mapping = aes(x = gamma_var_rate_comp, y = gamma_var_rate_agg), show.legend = F, size = 2, min.segment.length = 0.15, force_pull = 1.2, max.time = 1) +
  scale_x_log10(limits = c(0.001, NA), breaks = c(0.001, .01, .1)) + 
  scale_y_log10() +
  #geom_text_repel(size = 2.5) +
  scale_color_manual(values = pal, drop = FALSE) +
  theme(legend.position = "none") +
  coord_fixed() +
  labs(x = "Compositional variability",
       y = "Aggregate variability",
       color = "Organism group")
ggsave(filename = "figs/comp_agg_compare.png",plot = comp_agg_fig, bg = "white", width = 6, height = 6, dpi = 600)
 



s_rho <- cor(comp_agg_stab$phi_var_comp, comp_agg_stab$phi_var_agg, use = "pairwise.complete", method = "spearman")

phi_compare <- na.omit(comp_agg_stab) %>% 
  ggplot(aes(y = phi_var_agg,
             x = phi_var_comp,  label = lter_site, 
             color = organism_group, group = paste(lter_site, `organism_group`))) +
  geom_abline(slope = 1, intercept = 0, alpha = 0.25, linetype = "dashed") +
  geom_point(size = 3, alpha = 0.5) +
  scale_color_manual(values = pal, drop = FALSE) +
  geom_text_repel(show.legend = F, size = 3) +
  labs(color = "Organism group",
       y = expression(paste("Agg. Spatial Synchrony (",phi,")")),
       x = expression(paste("Comp. Spatial Synchrony (",BD[phi],")"))) +
  coord_fixed() +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8)) + 
  theme(legend.position = "right") +
  annotate("text", x = .75, y = 0.1, size = 5, label = expression(paste(rho, "= 0.53")))
ggsave("figs/phi_comparison.png",plot = phi_compare, bg = "white", width = 6, height = 3/4*6, dpi = 600)

phi_compare_fig <- comp_agg_fig + phi_compare + 
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", nrow = 1)
# ggsave("figs/agg_comp_panel.png", plot = phi_compare_fig, bg = "white", width = 8, height = 6, dpi = 600)

# 5. PARTITIONING DSRs --------------------------------------------------------
comp_mod <- psem(
  lm(gamma_var ~ phi_var + alpha_var + gamma_div_mean, data = metacom_divstab_comp_dat),
  lm(phi_var ~ beta_div_mean, data = metacom_divstab_comp_dat),
  lm(alpha_var ~ alpha_div_mean, data = metacom_divstab_comp_dat),
  data = metacom_divstab_comp_dat
)

summary(comp_mod)
plot(
  comp_mod,
  return = FALSE,
  node_attrs = data.frame(shape = "rectangle", color = "black", fillcolor = "white"),
  edge_attrs = data.frame(style = "solid", color = "black"),
  ns_dashed = T,
  alpha = 0.05,
  show = "std",
  digits = 3,
  add_edge_label_spaces = TRUE)

agg_mod <- psem(
  lm(gamma_var ~ phi_var + alpha_var + gamma_div_mean, data = metacom_divstab_agg_dat),
  lm(phi_var ~ beta_div_mean, data = metacom_divstab_agg_dat),
  lm(alpha_var ~ alpha_div_mean, data = metacom_divstab_agg_dat),
  data = metacom_divstab_agg_dat
)

summary(agg_mod)
plot(
  agg_mod,
  return = FALSE,
  node_attrs = data.frame(shape = "rectangle", color = "black", fillcolor = "white"),
  edge_attrs = data.frame(style = "solid", color = "black"),
  ns_dashed = T,
  alpha = 0.05,
  show = "std",
  digits = 3,
  add_edge_label_spaces = TRUE)

dsr_ag <- metacom_divstab_agg_dat %>% 
  rename(cv_gamma = gamma_var, 
         phi = phi_var,
         cv_alpha = alpha_var) %>% 
  select(dataset_id, cv_gamma, phi, cv_alpha, alpha_div_mean, beta_div_mean, gamma_div_mean)
dsr_com <- metacom_divstab_comp_dat %>% 
  rename(bd_gamma = gamma_var, 
         bd_phi = phi_var,
         bd_alpha = alpha_var) %>% 
  select(dataset_id, bd_gamma, bd_phi, bd_alpha, alpha_div_mean, beta_div_mean, gamma_div_mean)
dsr_tot <- left_join(dsr_ag, dsr_com)

write_csv(dsr_tot, file = "dsr_results_table.csv")


tot_mod <- psem(
  lm(cv_gamma ~ cv_alpha + phi + bd_alpha + bd_phi + bd_gamma + gamma_div_mean, data = dsr_tot),
  lm(bd_gamma ~ bd_alpha + bd_phi + gamma_div_mean, data = dsr_tot),
  lm(gamma_div_mean ~ alpha_div_mean + beta_div_mean, data = dsr_tot),
  lm(phi ~ beta_div_mean + bd_phi, data = dsr_tot),
  lm(cv_alpha ~ alpha_div_mean + bd_alpha, data = dsr_tot),
  lm(bd_phi ~ beta_div_mean, data = dsr_tot),
  lm(bd_alpha ~ alpha_div_mean, data = dsr_tot),
  data = dsr_tot
)
summary(tot_mod)
AIC_psem(tot_mod)
sem_coefs <- coefs(tot_mod)
rsquared(tot_mod)

sem_path <- plot(
  tot_mod,
  return = T,
  #node_attrs = data.frame(shape = "rectangle", color = "black", fillcolor = "white"),
  #edge_attrs = data.frame(style = "solid", color = "black"),
  ns_dashed = T,
  alpha = 0.05,
  show = "std",
  digits = 3,
  add_edge_label_spaces = TRUE)

sem_path$edges_df <- sem_path$edges_df %>% 
  mutate(color = ifelse(as.numeric(label) < 0, "red", "blue"),
         penwidth = abs(coefs(tot_mod)$Std.Estimate) * 10)

sem_path$edges_df

get_global_graph_attr_info(sem_path)
sem_path$global_attrs[7, 2] <- "false"
sem_path$global_attrs[1, 2] <- "dot"
?layout_nodes_w_string
sem_path$nodes_df$label <- c(
  "Agg. Metacommunity \nVariability (CV_gamma)",
  "Comp. Metacommunity \nVariability (BD_gamma)",
  "Avg. Gamma Diversity",
  "Agg. Spatial \nSynchrony (phi)",
  "Agg. Local \nVariability (CV_alpha)",
  "Comp. Spatial \nSynchrony (BD_phi)",
  "Comp. Local \nVariability (BD_alpha)",
  "Avg. Alpha Diversity",
  "Avg. Beta Diversity"
)
sem_path$edges_df
sem_path$global_attrs[5,2] <- 10
sem_path$global_attrs[14,2] <- 12

# remove non significant edges
sem_path <- sem_path %>%
  select_edges(conditions = style != "solid") %>%
  delete_edges_ws()

sem_plot <- DiagrammeR::render_graph(sem_path)
sem_plot

# indirect paths
0.991*0.089
0.594 * 0.089
-.62*.525

-.407*.683
.398*.835

# https://rpubs.com/tjmahr/sem_diagrammer
# 6. ENVIRONMENTS AND TRAITS -------------------------------------------------


