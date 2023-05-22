library(here)
library(broom)
library(broom.mixed)
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
library(rstanarm)
library(bayesplot)
library(tidybayes)

theme_set(theme_base() + 
            theme(plot.background = element_blank(),
                  legend.title = element_text(size = 10),
                  legend.text = element_text(size = 8)))

pal <- colorspace::darken(RColorBrewer::brewer.pal(n = 10, name = "Set3"), amount = .2)

set.seed(39759)

options(mc.cores = parallel::detectCores())

# 1. IMPORT DATA SETS ------------------------------------------------------------
#source("analysis_wrapper.R")
metacom_var <- read_csv(here("data/L4_metacommunity_variability_analysis_results_2023-04-19.csv"))
local_var <- read_csv(here("data/L4_local_variability_analysis_results_2023-04-19.csv"))
env_var <- read_csv(here("data/lter_centroid_satdata.csv"))
data_list <- read_csv(here("data/L3_DATA_list.csv")) %>% 
  filter(l3_filename != "L3-mcr-fish-castorani.csv") 
data_list$organism_group[which(data_list$organism_group=="coral")] <- "invertebrates"

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
(p_val <- as.character(round(local_div_stab_comp_alpha_fit$p.value,2))) # 0.16
(r2 <- as.character(round(local_div_stab_comp_alpha_fit$r.squared,2))) # 0.


local_div_stab_agg_alpha_mod <- lm(CV ~ alpha_div_scaled, data = local_dataset_for_mods)
local_div_stab_agg_alpha_fit <- glance(local_div_stab_agg_alpha_mod)
summary(local_div_stab_agg_alpha_mod)
(p_val <- as.character(round(local_div_stab_agg_alpha_fit$p.value,2))) # 0.
(r2 <- as.character(round(local_div_stab_agg_alpha_fit$r.squared,2))) # 0.09



# uncorrelated random intercept and slope
local_comp_mod_lmm.0 <- lmer(BD ~ alpha_div_scaled + (alpha_div_scaled||dataset_id), data = local_dataset_for_mods)
# correlated random intercept and slope
local_comp_mod_lmm.1 <- lmer(BD ~ alpha_div_scaled + (alpha_div_scaled|dataset_id), data = local_dataset_for_mods)
#random intercept fixed mean
local_comp_mod_lmm.2 <- lmer(BD ~ alpha_div_scaled + (1|dataset_id), data = local_dataset_for_mods)

AIC(local_comp_mod_lmm.0, local_comp_mod_lmm.1, local_comp_mod_lmm.2)
local_comp_mod_lmm <- local_comp_mod_lmm.1
summary(local_comp_mod_lmm)

plot(local_comp_mod_lmm)
coef(local_comp_mod_lmm)


# uncorrelated random intercept and slope
local_agg_mod_lmm.0 <- lmer(CV ~ alpha_div_scaled + (alpha_div_scaled||dataset_id), data = local_dataset_for_mods)
# correlated random intercept and slope
local_agg_mod_lmm.1 <- lmer(CV ~ alpha_div_scaled + (alpha_div_scaled|dataset_id), data = local_dataset_for_mods)
#random intercept fixed mean
local_agg_mod_lmm.2 <- lmer(CV ~ alpha_div_scaled + (1|dataset_id), data = local_dataset_for_mods)
AIC(local_agg_mod_lmm.0, local_agg_mod_lmm.1, local_agg_mod_lmm.2)
local_agg_mod_lmm <- local_agg_mod_lmm.1
summary(local_agg_mod_lmm)
plot(local_agg_mod_lmm)
coef(local_agg_mod_lmm)


# make re figs
ggsave(filename = "figs/local_agg_randomeffects.png", 
       plot = plot_model(local_agg_mod_lmm, type = "re"),
       width = 7, height = 6, dpi = 1000, bg = "white")
ggsave(filename = "figs/local_comp_randomeffects.png", 
       plot = plot_model(local_comp_mod_lmm, type = "re"),
       width = 7, height = 6, dpi = 1000, bg = "white")
tab_model(local_agg_mod_lmm, local_comp_mod_lmm)

(r2m_comp <- r.squaredGLMM(local_comp_mod_lmm)[1])
(r2c_comp <- r.squaredGLMM(local_comp_mod_lmm)[2])

(r2m_agg <- r.squaredGLMM(local_agg_mod_lmm)[1])
(r2c_agg <- r.squaredGLMM(local_agg_mod_lmm)[2])

### make figs

R2m_expression <- expression(paste(" ", R[m]^2 , "= ", 0.013))
R2c_expression <- expression(paste(" ", R[c]^2 , "= ", 0.727))

local_divstab_comp_fig <- local_dataset_for_mods %>% 
                          
  ggplot(aes(x = alpha_div_scaled, y = BD)) + 
  geom_point(mapping = aes(group = dataset_id, color = organism_group, shape = organism_group), alpha = 0.5) + 
  geom_smooth(mapping = aes(group= dataset_id, color = organism_group), method = "lm", linewidth=0.5, se = F, show.legend = FALSE) + 
  geom_smooth(method = "lm", se = F, color = "black", linewidth = 1.5) +
  labs(x = expression(paste("Mean ", alpha, "-diversity (z-score)")), y = expression(paste("Comp. ", alpha, "-variability (", BD[alpha]^h,")")), color = "Organism group", shape = "Organism group") + 
  scale_color_manual(values = pal) +
  #scale_y_log10() +
  annotate("text", x = 2.5, y = 0.8, label = R2m_expression) +
  annotate("text", x = 2.5, y = 0.7, label = R2c_expression)

local_divstab_comp_fig


R2m_expression <- expression(paste(" ", R[m]^2 , "= ", 0.022))
R2c_expression <- expression(paste(" ", R[c]^2 , "= ", 0.744))

local_divstab_agg_fig <- local_dataset_for_mods %>% 
  
  ggplot(aes(x = alpha_div_scaled, y = CV)) + 
  geom_point(mapping = aes(group = dataset_id, color = organism_group, shape = organism_group), alpha = 0.5) + 
  geom_smooth(mapping = aes(group = dataset_id, color = organism_group), method = "lm", size =0.5, se = F, show.legend = FALSE) + 
  geom_smooth(method = "lm", se = F, color = "black", size = 1.5) +
  labs(x = expression(paste("Mean ", alpha, "-diversity (z-score)")), y = expression(paste("Agg. ", alpha, "-variability (CV)")), color = "Organism group", shape = "Organism group") + 
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

# div_stab_comp_gamma_mod <- (lm(gamma_var ~ gamma_div_mean, data = metacom_divstab_comp_dat))
# div_stab_comp_beta_mod <- (lm(phi_var ~ beta_div_mean, data = metacom_divstab_comp_dat))
# div_stab_comp_alpha_mod <- (lm(alpha_var ~ alpha_div_mean, data = metacom_divstab_comp_dat))
# div_stab_comp_alpha_gamma_mod <- (lm(gamma_var ~ alpha_div_mean, data = metacom_divstab_comp_dat))
# div_stab_comp_beta_gamma_mod <- (lm(gamma_var ~ beta_div_mean, data = metacom_divstab_comp_dat))
# 
# summary(div_stab_comp_alpha_mod)
# summary(div_stab_comp_beta_mod)
# summary(div_stab_comp_gamma_mod)
# summary(div_stab_comp_alpha_gamma_mod)
# summary(div_stab_comp_beta_gamma_mod)
# 
# dsrc_g  <- (lmer(gamma_var ~ gamma_div_mean + (gamma_div_mean|lter_site), data = metacom_divstab_comp_dat))
# dsrc_ga <- (lmer(gamma_var ~ alpha_div_mean + (alpha_div_mean|lter_site), data = metacom_divstab_comp_dat))
# dsrc_gb <- (lmer(gamma_var ~ beta_div_mean + (beta_div_mean|lter_site), data = metacom_divstab_comp_dat))
# dsrc_b  <- (lmer(phi_var ~ beta_div_mean + (beta_div_mean|lter_site), data = metacom_divstab_comp_dat))
# dsrc_a  <- (lmer(alpha_var ~ alpha_div_mean + (alpha_div_mean|lter_site), data = metacom_divstab_comp_dat))



dsrc_g_b <- stan_lmer(gamma_var ~ gamma_div_mean + (gamma_div_mean|lter_site), data = metacom_divstab_comp_dat, iter = 5000, adapt_delta = 0.99)
dsrc_ga_b <- stan_lmer(gamma_var ~ alpha_div_mean + (alpha_div_mean|lter_site), data = metacom_divstab_comp_dat, iter = 5000, adapt_delta = 0.99)
dsrc_gb_b <- stan_lmer(gamma_var ~ beta_div_mean + (beta_div_mean|lter_site), data = metacom_divstab_comp_dat, iter = 5000, adapt_delta = 0.99)
dsrc_b_b  <- stan_lmer(phi_var ~ beta_div_mean + (beta_div_mean|lter_site), data = metacom_divstab_comp_dat, iter = 5000, adapt_delta = 0.99)
dsrc_a_b  <- stan_lmer(alpha_var ~ alpha_div_mean + (alpha_div_mean|lter_site), data = metacom_divstab_comp_dat, iter = 5000, adapt_delta = 0.99)

# div_stab_agg_gamma_mod <- (lm(gamma_var ~ gamma_div_mean, data = metacom_divstab_agg_dat))
# div_stab_agg_beta_mod <- (lm(phi_var ~ beta_div_mean, data = metacom_divstab_agg_dat))
# div_stab_agg_alpha_mod <- (lm(alpha_var ~ alpha_div_mean, data = metacom_divstab_agg_dat))
# div_stab_agg_alpha_gamma_mod <- (lm(gamma_var ~ alpha_div_mean, data = metacom_divstab_agg_dat))
# div_stab_agg_beta_gamma_mod <- (lm(gamma_var ~ beta_div_mean, data = metacom_divstab_agg_dat))
# 
# summary(div_stab_agg_alpha_mod)
# summary(div_stab_agg_beta_mod)
# summary(div_stab_agg_gamma_mod)
# summary(div_stab_agg_alpha_gamma_mod)
# summary(div_stab_agg_beta_gamma_mod)
# 
# dsra_g  <- (lmer(gamma_var ~ gamma_div_mean + (gamma_div_mean|lter_site), data = metacom_divstab_agg_dat))
# dsra_ga <- (lmer(gamma_var ~ alpha_div_mean + (alpha_div_mean|lter_site), data = metacom_divstab_agg_dat))
# dsra_gb <- (lmer(gamma_var ~ beta_div_mean + (beta_div_mean|lter_site), data = metacom_divstab_agg_dat))
# dsra_b  <- (lmer(phi_var ~ beta_div_mean + (beta_div_mean|lter_site), data = metacom_divstab_agg_dat))
# dsra_a  <- (lmer(alpha_var ~ alpha_div_mean + (alpha_div_mean|lter_site), data = metacom_divstab_agg_dat))

dsra_g_b  <- stan_lmer(gamma_var ~ gamma_div_mean + (gamma_div_mean|lter_site), data = metacom_divstab_agg_dat, iter = 5000, adapt_delta = 0.99)
dsra_ga_b <- stan_lmer(gamma_var ~ alpha_div_mean + (alpha_div_mean|lter_site), data = metacom_divstab_agg_dat, iter = 5000, adapt_delta = 0.99)
dsra_gb_b <- stan_lmer(gamma_var ~ beta_div_mean + (beta_div_mean|lter_site), data = metacom_divstab_agg_dat, iter = 5000, adapt_delta = 0.99)
dsra_b_b  <- stan_lmer(phi_var ~ beta_div_mean + (beta_div_mean|lter_site), data = metacom_divstab_agg_dat, iter = 5000, adapt_delta = 0.99)
dsra_a_b  <- stan_lmer(alpha_var ~ alpha_div_mean + (alpha_div_mean|lter_site), data = metacom_divstab_agg_dat, iter = 5000, adapt_delta = 0.99)

median(bayes_R2(dsrc_g_b ))
median(bayes_R2(dsrc_ga_b))
median(bayes_R2(dsrc_gb_b))
median(bayes_R2(dsrc_b_b ))
median(bayes_R2(dsrc_a_b ))
median(bayes_R2(dsra_g_b ))
median(bayes_R2(dsra_ga_b))
median(bayes_R2(dsra_gb_b))
median(bayes_R2(dsra_b_b ))
median(bayes_R2(dsra_a_b ))

prior_summary(dsrc_g_b)
prior_summary(dsrc_ga_b)
prior_summary(dsrc_gb_b)
prior_summary(dsrc_b_b )
prior_summary(dsrc_a_b )
prior_summary(dsra_g_b )
prior_summary(dsra_ga_b)
prior_summary(dsra_gb_b)
prior_summary(dsra_b_b )
prior_summary(dsra_a_b )

pdf(file = "figs/posterior_fits.pdf", width = 12, height = 6)
par(mfrow = c(2, 5))
boxplot(posterior_predict(dsrc_g_b), pch = NA, ylab = "Comp. Gamma Var")
points(metacom_divstab_comp_dat$gamma_var, pch = 16, col = "red")
boxplot(posterior_predict(dsrc_gb_b), pch = NA, ylab = "Comp. Gamma Var")
points(metacom_divstab_comp_dat$gamma_var, pch = 16, col = "red")
boxplot(posterior_predict(dsrc_ga_b), pch = NA, ylab = "Comp. Gamma Var")
points(metacom_divstab_comp_dat$gamma_var, pch = 16, col = "red")
boxplot(posterior_predict(dsrc_b_b), pch = NA, ylab = "Comp. Phi")
points(metacom_divstab_comp_dat$phi_var, pch = 16, col = "red")
boxplot(posterior_predict(dsrc_a_b), pch = NA, ylab = "Comp. Alpha Var")
points(metacom_divstab_comp_dat$alpha_var, pch = 16, col = "red")

boxplot(posterior_predict(dsra_g_b), pch = NA, ylab = "Agg. Gamma Var")
points(metacom_divstab_agg_dat$gamma_var, pch = 16, col = "red")
boxplot(posterior_predict(dsra_gb_b), pch = NA, ylab = "Agg. Gamma Var")
points(metacom_divstab_agg_dat$gamma_var, pch = 16, col = "red")
boxplot(posterior_predict(dsra_ga_b), pch = NA, ylab = "Agg. Gamma Var")
points(metacom_divstab_agg_dat$gamma_var, pch = 16, col = "red")
boxplot(posterior_predict(dsra_b_b), pch = NA, ylab = "Agg. Phi")
points(metacom_divstab_agg_dat$phi_var, pch = 16, col = "red")
boxplot(posterior_predict(dsra_a_b), pch = NA, ylab = "Agg. Alpha Var")
points(metacom_divstab_agg_dat$alpha_var, pch = 16, col = "red")
dev.off()

prior_summary(dsra_g_b)




plot_model(dsra_gb_b, type = "re")
plot_model(dsrc_g_b, type = "re")

b0 <- tidy(dsra_gb_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsra_gb_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]

metacom_divstab_agg_dat %>% 
  add_epred_draws(dsra_gb_b, ndraws = 200, re_formla = NA) %>% 
  ggplot(aes(x = beta_div_mean, y = gamma_var)) + 
  geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
  geom_abline(intercept = b0, slope = b1, size = 1, color = "blue") +
  geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
  #geom_label_repel(size = 2) +
  labs(x = expression(paste("Mean ", gamma, "-diversity")),
       y = expression(paste("Comp. ", gamma, "-variability")), # (", BD[gamma]^h, ") for metric
       color = "Organism group", shape = "Organism group")  +
  #scale_x_log10() +
  scale_color_manual(values = pal, drop = FALSE)

tab_model(dsrc_g_b, dsrc_gb_b, dsrc_ga_b,  dsrc_b_b, dsrc_a_b,  
          pred.labels = c("Intercept", "Mean Gamma-Diversity", 
                          "Mean Beta-Diversity", "Mean Alpha-Diversity"), 
          dv.labels = c("Regional (BD_gamma)", "Regional (BD_gamma)",
                        "Regional (BD_gamma)","Spatial Synchrony (BD_phi)", 
                        "Local (BD_alpha)"),
          file = "tables/DSRs_comp_models.html")
tab_model(dsra_g_b, dsra_gb_b, dsra_ga_b,  dsra_b_b, dsra_a_b, show.re.var = T, 
          pred.labels = c("Intercept", "Mean Gamma-Diversity", 
                          "Mean Beta-Diversity", "Mean Alpha-Diversity"), 
          dv.labels = c("Regional (CV_gamma)", "Regional (CV_gamma)",
                        "Regional (CV_gamma)","Spatial Synchrony (phi)", 
                        "Local (CV_alpha)"), 
          file = "tables/DSRs_agg_models.html")

comp_var_plot <- plot_models(dsrc_g_b, dsrc_gb_b, dsrc_ga_b, dsrc_b_b, dsrc_a_b, 
            title = "Compositional Variability", 
            legend.title = "Variability Partition", 
            m.labels = c("Regional (BD_gamma)", "Regional (BD_gamma)",
                         "Regional (BD_gamma)","Spatial Synchrony (BD_phi)", 
                         "Local (BD_alpha)"), show.values = TRUE, 
            axis.labels = c("Mean Alpha-Diversity", "Mean Beta-Diversity",
                            "Mean Gamma-Diversity"), show.p = F)
agg_var_plot <- plot_models(dsra_g_b, dsra_gb_b, dsra_ga_b, dsra_b_b, dsra_a_b,
            title = "Aggregate Variability", 
            legend.title = "Variability Partition", 
            m.labels = c("Regional (CV_gamma)", "Regional (CV_gamma)",
                         "Regional (CV_gamma)","Spatial Synchrony (phi)", 
                         "Local (CV_alpha)"), show.values = TRUE, 
            axis.labels = c("Mean Alpha-Diversity", "Mean Beta-Diversity", 
                            "Mean Gamma-Diversity"), show.p = F)
var_estimates_plot <- agg_var_plot + comp_var_plot + patchwork::plot_layout(ncol = 1)
var_estimates_plot
ggsave(filename = "figs/dsr_estimates.png", width = 8, height = 8, dpi = 1000, bg = "white")


# regional compositional variability

# div_stab_comp_gamma_fit <- glance(div_stab_comp_gamma_mod)
# (p_val <- as.character(round(div_stab_comp_gamma_fit$p.value,2))) # 0.05
# (r2 <- as.character(round(div_stab_comp_gamma_fit$r.squared,2))) # 0.17
# p_expression <- expression(paste(" ", p , "= ", 0.11))
# R2_expression <- expression(paste(" ", R^2 , "= ", 0.09))
b0 <- tidy(dsrc_g_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsrc_g_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]

(div_stab_gamma_h <- metacom_divstab_comp_dat %>% 
    add_epred_draws(dsrc_g_b, ndraws = 200, re_formla = NA) %>% 
    ggplot(aes(x = gamma_div_mean, y = gamma_var, label = lter_site)) +
    geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
    #geom_abline(intercept = b0, slope = b1, size = 1, color = "black") +
    #stat_smooth(method = "lm", se = T, size = 1, color = "black") +
    geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
    #geom_label_repel(size = 2) +
    labs(x = expression(paste("Mean ", gamma, "-diversity")),
         y = expression(paste("Comp. ", gamma, "-variability")), # (", BD[gamma]^h, ") for metric
         color = "Organism group", shape = "Organism group")  +
    #scale_x_log10() +
    scale_color_manual(values = pal, drop = FALSE) 
    #annotate("text", x = 12, y = .6, label = expression(paste(" ", p , "= ", 0.11))) +
    #annotate("text", x = 12, y = .55, label = expression(paste(" ", R^2 , "= ", 0.09)))
  
  #ggsave("ESA_2019/figs/variability_alpha-gamma.png", width = 6, height = 4, units = "in", dpi = 600)
)

# div_stab_comp_alpha_fit <- glance(div_stab_comp_alpha_mod)
# (p_val <- as.character(round(div_stab_comp_alpha_fit$p.value,2))) # 0.86
# (r2 <- as.character(round(div_stab_comp_alpha_fit$r.squared,2))) # 0.
# p_expression <- expression(paste(" ", p , "= ", 0.22))
# R2_expression <- expression(paste(" ", R^2 , "= ", 0.05))

b0 <- tidy(dsrc_a_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsrc_a_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]

(div_stab_alpha_h <- metacom_divstab_comp_dat %>% 
    add_epred_draws(dsrc_a_b, ndraws = 200, re_formla = NA) %>% 
    ggplot(aes(x = alpha_div_mean, y = alpha_var, label = lter_site)) +
    geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
    #geom_abline(intercept = b0, slope = b1, size = 1, color = "black") +
    # stat_smooth(method = "lm", se = T, size = 1, color = "black") +
    geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
    #geom_label_repel(size = 2) +
    labs(x = expression(paste("Mean ", alpha, "-diversity")),
         y = expression(paste("Comp. ", alpha, "-variability")),
         color = "Organism group", shape = "Organism group") +
    scale_color_manual(values = pal, drop = FALSE)
    #scale_x_log10() +
    # annotate("text", x = 25, y = .2, label = expression(paste(" ", p , "= ", 0.22))) +
    # annotate("text", x = 25, y = .15, label = expression(paste(" ", R^2 , "= ", 0.05)))
 
)


# div_stab_alpha_gamma_fit <- glance(div_stab_comp_alpha_gamma_mod)
# div_stab_alpha_gamma_fit$p.value 
# div_stab_alpha_gamma_fit$r.squared 
# p_expression <- expression(paste(" ", p , "= ", 0.011))
# R2_expression <- expression(paste(" ", R^2 , "= ", 0.217))

b0 <- tidy(dsrc_ga_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsrc_ga_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]

(div_stab_alpha_gamma_h <- metacom_divstab_comp_dat %>% 
    add_epred_draws(dsrc_ga_b, ndraws = 200, re_formula = NA) %>% 
    ggplot(aes(x = alpha_div_mean, y = gamma_var, label = lter_site)) +
    geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
    geom_abline(intercept = b0, slope = b1, size = 1, color = "black") +
    # stat_smooth(method = "lm", se = T, size = 1, color = "black") +
    geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
    #geom_label_repel(size = 2) +
    labs(x = expression(paste("Mean ", alpha, "-diversity")),
         y = expression(paste("Comp. ", gamma, "-variability")),
         color = "Organism group", shape = "Organism group") +
    scale_color_manual(values = pal, drop = FALSE)
    #scale_x_log10()  +
    # annotate("text", x = 25, y = .1, label = expression(paste(" ", p , "= ", 0.011))) +
    # annotate("text", x = 25, y = .05, label = expression(paste(" ", R^2 , "= ", 0.217)))
)

# div_stab_beta_gamma_fit <- glance(div_stab_comp_beta_gamma_mod)
# div_stab_beta_gamma_fit$p.value # 0.2
# div_stab_beta_gamma_fit$r.squared
# p_expression <- expression(paste(" ", p , "= ", 0.42))
# R2_expression <- expression(paste(" ", R^2 , "= ", 0.02))

b0 <- tidy(dsrc_gb_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsrc_gb_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]

(div_stab_beta_h <- metacom_divstab_comp_dat %>% 
    add_epred_draws(dsrc_gb_b, ndraws = 200, re_formula = NA) %>% 
    ggplot(aes(x = beta_div_mean, y = gamma_var, label = lter_site)) +
    geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
    #geom_abline(intercept = b0, slope = b1, size = 1, color = "black") +
    #stat_smooth(method = "lm", se = T, size = 1, color = "black") +
    geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
    scale_color_manual(values = pal, drop = FALSE) +
    #geom_label_repel(size = 2) +
    labs(x = expression(paste("Mean ", beta, "-diversity")),
         y = expression(paste("Comp. ", gamma, "-variability")),
         color = "Organism group", shape = "Organism group")
    # annotate("text", x = 1.75, y = .5, label = expression(paste(" ", p , "= ", 0.42))) +
    # annotate("text", x = 1.75, y = .45, label = expression(paste(" ", R^2 , "= ", 0.02)))
)

# regional aggregate variability
# div_stab_gamma_agg_fit <- glance(div_stab_agg_gamma_mod)
# div_stab_gamma_agg_fit$p.value # 0.202
# div_stab_gamma_agg_fit$r.squared # 0.08
# p_expression <- expression(paste(" ", p , "= ", 0.25))
# R2_expression <- expression(paste(" ", R^2 , "= ", 0.049))

b0 <- tidy(dsra_g_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsra_g_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]

(div_stab_gamma_agg <- metacom_divstab_agg_dat %>% 
    add_epred_draws(dsra_g_b, ndraws = 200, re_formula = NA) %>% 
    ggplot(aes(x = gamma_div_mean, y = gamma_var, label = lter_site)) +
    geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
    #geom_abline(intercept = b0, slope = b1, size = 1, color = "black") +
    #stat_smooth(method = "lm", se = T, size = 1, color = "black") +
    geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
    scale_color_manual(values = pal, drop = FALSE) +
    #geom_label_repel(size = 2) +
    labs(x = expression(paste("Mean ", gamma, "-diversity")),
         y = expression(paste("Agg. ", gamma, "-variability")),
         color = "Organism group", shape = "Organism group") 
    # annotate("text", x = 60, y = .75, label = expression(paste(" ", p , "= ", 0.25))) +
    # annotate("text", x = 60, y = .65, label = expression(paste(" ", R^2 , "= ", 0.049)))
)

# div_stab_alpha_agg_fit <- glance(div_stab_agg_alpha_mod)
# div_stab_alpha_agg_fit$p.value # 0.17
# div_stab_alpha_agg_fit$r.squared # 0.0927
# p_expression <- expression(paste(" ", p , "= ", 0.45))
# R2_expression <- expression(paste(" ", R^2 , "= ", 0.022))
b0 <- tidy(dsra_a_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsra_a_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]

(div_stab_alpha_agg <- metacom_divstab_agg_dat %>% 
    add_epred_draws(dsra_a_b, ndraws = 200, re_formula = NA) %>% 
    ggplot(aes(x = alpha_div_mean, y = alpha_var, label = lter_site)) +
    geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
    #geom_abline(intercept = b0, slope = b1, size = 1, color = "black") +
    #stat_smooth(method = "lm", se = T, size = 1, color = "black") +
    geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
    scale_color_manual(values = pal, drop = FALSE) +
    #geom_label_repel(size = 2) +
    labs(x = expression(paste("Mean ", alpha, "-diversity")),
         y = expression(paste("Agg. ", alpha, "-variability")),
         color = "Organism group", shape = "Organism group")
    # scale_x_log10() +
    # annotate("text", x = 25, y = .9, label = expression(paste(" ", p , "= ", 0.45))) +
    # annotate("text", x = 25, y = .85, label = expression(paste(" ", R^2 , "= ", 0.022)))
)

# div_stab_alpha_gamma_agg_fit <- glance(div_stab_agg_alpha_gamma_mod)
# div_stab_alpha_gamma_agg_fit$p.value # 0.927
# div_stab_alpha_gamma_agg_fit$r.squared # 0
# p_expression <- expression(paste(" ", p , "= ", 0.65))
# R2_expression <- expression(paste(" ", R^2 , "= ", 0.008))
b0 <- tidy(dsra_ga_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsra_ga_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]

(div_stab_alpha_gamma_agg <- metacom_divstab_agg_dat %>% 
    add_epred_draws(dsra_ga_b, ndraws = 200, re_formula = NA) %>% 
    ggplot(aes(x = alpha_div_mean, y = gamma_var, label = lter_site)) +
    #stat_smooth(method = "lm", se = T, size = 1, color = "black") +
    geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
    #geom_abline(intercept = b0, slope = b1, size = 1, color = "black") +
    geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
    scale_color_manual(values = pal, drop = FALSE) +
    #geom_label_repel(size = 2) +
    labs(x = expression(paste("Mean ", alpha, "-diversity")),
         y = expression(paste("Agg. ", gamma, "-variability")),
         color = "Organism group", shape = "Organism group")
    #scale_x_log10() +
    # annotate("text", x = 25, y = .8, label = expression(paste(" ", p , "= ", 0.65))) +
    # annotate("text", x = 25, y = .7, label = expression(paste(" ", R^2 , "= ", 0.008)))
)

# div_stab_beta_gamma_agg_fit <- glance(div_stab_agg_beta_gamma_mod)
# div_stab_beta_gamma_agg_fit$p.value # 0.0095
# div_stab_beta_gamma_agg_fit$r.squared # 0.29
# p_expression <- expression(paste(" ", p , "= ", 0.026))
# R2_expression <- expression(paste(" ", R^2 , "= ", 0.171))

b0 <- tidy(dsra_gb_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsra_gb_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]


(div_stab_beta_agg <- metacom_divstab_agg_dat %>% 
    add_epred_draws(dsra_ga_b, ndraws = 200, re_formula = NA) %>% 
    ggplot(aes(x = beta_div_mean, y = gamma_var)) +
    geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
    geom_abline(intercept = b0, slope = b1, size = 1, color = "black") +
    # stat_smooth(method = "lm", se = T, size = 1, color = "black") +
    geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
    scale_color_manual(values = pal, drop = FALSE) +
    #geom_label_repel(size = 2) +
    labs(x = expression(paste("Mean ", beta, "-diversity")),
         y = expression(paste("Agg. ", gamma, "-variability")),
         color = "Organism group", shape = "Organism group")
    # annotate("text", x = 4.5, y = .7, label = expression(paste(" ", p , "= ", 0.026))) +
    # annotate("text", x = 4.5, y = .6, label = expression(paste(" ", R^2 , "= ", 0.171)))
)

# combine figs
div_stab_multipanel_fig <- div_stab_gamma_agg + div_stab_gamma_h +
  div_stab_beta_agg + div_stab_beta_h +
  div_stab_alpha_gamma_agg + div_stab_alpha_gamma_h +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect", nrow = 3)
ggsave(filename = "figs/diversity_variability.png", plot = div_stab_multipanel_fig, width = 2.5*4, height = 3*3, units = "in", dpi = 1000, bg = "white")



# 4. COMPARE BETA DIVERSITY WITH PHI, SYNCHRONY
# summary(div_stab_comp_beta_mod)
# div_stab_phi_fit <- glance(div_stab_comp_beta_mod)
# div_stab_phi_fit$p.value # 2 e-4
# div_stab_phi_fit$r.squared # 0.49

b0 <- tidy(dsrc_b_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsrc_b_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]

(div_stab_phi_h <- metacom_divstab_comp_dat %>% 
    add_epred_draws(dsrc_b_b, ndraws = 200, re_formula = NA) %>% 
    ggplot(aes(x = beta_div_mean, y = phi_var)) +
    geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
    geom_abline(intercept = b0, slope = b1, size = 1, color = "black") +
    # stat_smooth(method = "lm", se = T, size = 1, color = "black") +
    geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
    scale_y_continuous(limits = c(0,1)) +
    scale_color_manual(values = pal, drop = FALSE) +
    #geom_label_repel(size = 2) +
    labs(x = expression(paste("Mean ", beta, "-diversity")),
         y = expression(paste("Comp. Synchrony (", BD[phi]^h, ")")),
         color = "Organism group", shape = "Organism group")
    # annotate("text", x = 4.75, y = 0.9, label = bquote(atop(paste(R^2, "= 0.38"),
    #                                                         "p = 0.0003")))
  #ggsave("ESA_2019/figs/variability_alpha-gamma.png", width = 6, height = 4, units = "in", dpi = 600)
)

# summary(div_stab_agg_beta_mod)
# div_stab_phi_agg_fit <- glance(div_stab_agg_beta_mod)
# div_stab_phi_agg_fit$p.value # 0.006
# div_stab_phi_agg_fit$r.squared # 0.32

b0 <- tidy(dsra_b_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsra_b_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]

(div_stab_phi_agg <- metacom_divstab_agg_dat %>% 
    add_epred_draws(dsra_b_b, ndraws = 200, re_formula = NA) %>% 
    ggplot(aes(x = beta_div_mean, y = phi_var)) +
    geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
    geom_abline(intercept = b0, slope = b1, size = 1, color = "black") +
    # stat_smooth(method = "lm", se = T, size = 1, color = "black") +
    geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
    scale_y_continuous(limits = c(0,1))+
    scale_color_manual(values = pal, drop = FALSE) +
    #geom_label_repel(size = 2) +
    labs(x = expression(paste("Mean ", beta, "-diversity")),
         y = expression(paste("Agg. Synchrony (",phi,")")),
         color = "Organism group", shape = "Organism group") 
    # annotate("text", x = 4.75, y = 0.9, label = bquote(atop(paste(R^2, "= 0.36"),
    #                                                         "p = 0.0006")))
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
             color = organism_group, shape = organism_group, group = paste(lter_site, `organism_group`))) +
  geom_abline(slope = 1, intercept = 0, alpha = 0.25, linetype = "dashed") +
  geom_point(size = 3) +
  scale_color_manual(values = pal, drop = FALSE) +
  geom_text_repel(show.legend = F, size = 3) +
  labs(color = "Organism group", shape = "Organism group",
       y = expression(paste("Agg. Spatial Synchrony (",phi,")")),
       x = expression(paste("Comp. Spatial Synchrony (",BD[phi],")"))) +
  coord_fixed() +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8)) + 
  scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8)) + 
  theme(legend.position = "right") +
  annotate("text", x = .75, y = 0.1, size = 5, label = expression(paste(rho, "= 0.53")))
ggsave("figs/phi_comparison.png",plot = phi_compare, bg = "white", width = 6, height = 3/4*6, dpi = 600)

phi_compare_fig <- comp_agg_fig + phi_compare + 
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", nrow = 1)
# ggsave("figs/agg_comp_panel.png", plot = phi_compare_fig, bg = "white", width = 8, height = 6, dpi = 600)

# 5. PARTITIONING DSRs --------------------------------------------------------

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


# https://rpubs.com/tjmahr/sem_diagrammer
# 6. ENVIRONMENTS AND TRAITS -------------------------------------------------


values_by_taxa_fig <- metacom_var %>% 
  filter(metric %in% c("gamma_var", "alpha_var", "phi_var")) %>% 
  mutate(metric = ifelse(metric == "gamma_var", "Metacommunity",
                         ifelse(metric == "phi_var", "Spatial synchrony", "Local"))) %>% 
  mutate(metric = factor(metric, levels = c("Metacommunity", "Spatial synchrony", "Local"))) %>% 
  mutate(variability_type = ifelse(variability_type == "agg", 
                                   "Aggregate", "Compositional")) %>% 
  ggplot(aes(x = organism_group, y = metric_value, fill = organism_group)) + 
  geom_boxplot() + 
  scale_fill_manual(values = pal, drop = FALSE) +
  facet_grid(metric ~ variability_type) +
  coord_flip() +
  theme(legend.position = "none") + 
  theme(panel.grid.major.x = element_line(color = "grey90")) +
  labs(x = "", y = "", fill = "Organism group")

ggsave(plot = values_by_taxa_fig, 
       filename = "figs/variability_taxa.png", 
       width = 6, height = 6, dpi = 1000, bg = "white")
