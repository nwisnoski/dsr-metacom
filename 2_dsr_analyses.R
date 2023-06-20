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
metacom_var <- read_csv(here("results/L4_metacommunity_variability_analysis_results_2023-06-20.csv"))
local_var <- read_csv(here("results/L4_local_variability_analysis_results_2023-06-20.csv"))
env_var <- read_csv(here("data/lter_centroid_satdata.csv"))
data_list <- read_csv(here("data/L3_DATA_list.csv"))
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

unique(local_var$dataset_id) 
unique(metacom_var$dataset_id) 

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

# uncorrelated random intercept and slope
local_comp_mod_lmm.0 <- lmer(BD ~ alpha_div_scaled + (alpha_div_scaled||dataset_id), data = local_dataset_for_mods)
# correlated random intercept and slope
local_comp_mod_lmm.1 <- lmer(BD ~ alpha_div_scaled + (alpha_div_scaled|dataset_id), data = local_dataset_for_mods)
# random intercept fixed mean
local_comp_mod_lmm.2 <- lmer(BD ~ alpha_div_scaled + (1|dataset_id), data = local_dataset_for_mods)
# nonlinear
local_comp_mod_lmm.3 <- lmer(BD ~ alpha_div_scaled + I(alpha_div_scaled^2) + (1|dataset_id), data = local_dataset_for_mods)
local_comp_mod_lmm.4 <- lmer(BD ~ alpha_div_scaled + I(alpha_div_scaled^2) + (alpha_div_scaled|dataset_id), data = local_dataset_for_mods)


AIC(local_comp_mod_lmm.0, local_comp_mod_lmm.1, local_comp_mod_lmm.2, local_comp_mod_lmm.3, local_comp_mod_lmm.4)
local_comp_mod_lmm <- local_comp_mod_lmm.1
summary(local_comp_mod_lmm)

plot(local_comp_mod_lmm)
coef(local_comp_mod_lmm)


# uncorrelated random intercept and slope
local_agg_mod_lmm.0 <- lmer(CV ~ alpha_div_scaled + (alpha_div_scaled||dataset_id), data = local_dataset_for_mods)
# correlated random intercept and slope
local_agg_mod_lmm.1 <- lmer(CV ~ alpha_div_scaled + (alpha_div_scaled|dataset_id), data = local_dataset_for_mods)
# random intercept fixed mean
local_agg_mod_lmm.2 <- lmer(CV ~ alpha_div_scaled + (1|dataset_id), data = local_dataset_for_mods)
# nonlinear
local_agg_mod_lmm.3 <- lmer(CV ~ alpha_div_scaled + I(alpha_div_scaled^2) + (1|dataset_id), data = local_dataset_for_mods)
local_agg_mod_lmm.4 <- lmer(CV ~ alpha_div_scaled + I(alpha_div_scaled^2) + (alpha_div_scaled|dataset_id), data = local_dataset_for_mods)

AIC(local_agg_mod_lmm.0, local_agg_mod_lmm.1, local_agg_mod_lmm.2, local_agg_mod_lmm.3, local_agg_mod_lmm.4)
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
tab_model(local_agg_mod_lmm, local_comp_mod_lmm, file = "tables/DSR_local_models.html")

(r2m_comp <- r.squaredGLMM(local_comp_mod_lmm)[1])
(r2c_comp <- r.squaredGLMM(local_comp_mod_lmm)[2])

(r2m_agg <- r.squaredGLMM(local_agg_mod_lmm)[1])
(r2c_agg <- r.squaredGLMM(local_agg_mod_lmm)[2])

### make figs

R2m_expression <- expression(paste(" ", R[m]^2 , "= ", 0.013))
R2c_expression <- expression(paste(" ", R[c]^2 , "= ", 0.726))

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


R2m_expression <- expression(paste(" ", R[m]^2 , "= ", 0.025))
R2c_expression <- expression(paste(" ", R[c]^2 , "= ", 0.715))

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

# bayesian mixed effects models
# compositional
dsrc_g_b <- stan_lmer(gamma_var ~ gamma_div_mean + (gamma_div_mean|lter_site), data = metacom_divstab_comp_dat, iter = 5000, adapt_delta = 0.99)
dsrc_ga_b <- stan_lmer(gamma_var ~ alpha_div_mean + (alpha_div_mean|lter_site), data = metacom_divstab_comp_dat, iter = 5000, adapt_delta = 0.99)
dsrc_gb_b <- stan_lmer(gamma_var ~ beta_div_mean + (beta_div_mean|lter_site), data = metacom_divstab_comp_dat, iter = 5000, adapt_delta = 0.99)
dsrc_b_b  <- stan_lmer(phi_var ~ beta_div_mean + (beta_div_mean|lter_site), data = metacom_divstab_comp_dat, iter = 5000, adapt_delta = 0.99)
dsrc_a_b  <- stan_lmer(alpha_var ~ alpha_div_mean + (alpha_div_mean|lter_site), data = metacom_divstab_comp_dat, iter = 5000, adapt_delta = 0.999)

# # nonlinear
# dsrc_g_b_2 <- stan_lmer(gamma_var ~ gamma_div_mean + I(gamma_div_mean^2) + (gamma_div_mean|lter_site), data = metacom_divstab_comp_dat, iter = 5000, adapt_delta = 0.99)
# dsrc_ga_b_2 <- stan_lmer(gamma_var ~ alpha_div_mean + I(alpha_div_mean^2) + (alpha_div_mean|lter_site), data = metacom_divstab_comp_dat, iter = 5000, adapt_delta = 0.99)
# dsrc_gb_b_2 <- stan_lmer(gamma_var ~ beta_div_mean + I(beta_div_mean^2) + (beta_div_mean|lter_site), data = metacom_divstab_comp_dat, iter = 5000, adapt_delta = 0.99)
# dsrc_b_b_2  <- stan_lmer(phi_var ~ beta_div_mean + I(beta_div_mean^2) + (beta_div_mean|lter_site), data = metacom_divstab_comp_dat, iter = 5000, adapt_delta = 0.99)
# dsrc_a_b_2  <- stan_lmer(alpha_var ~ alpha_div_mean + I(alpha_div_mean^2) + (alpha_div_mean|lter_site), data = metacom_divstab_comp_dat, iter = 5000, adapt_delta = 0.99)

# aggregate
dsra_g_b  <- stan_lmer(gamma_var ~ gamma_div_mean + (gamma_div_mean|lter_site), data = metacom_divstab_agg_dat, iter = 5000, adapt_delta = 0.99)
dsra_ga_b <- stan_lmer(gamma_var ~ alpha_div_mean + (alpha_div_mean|lter_site), data = metacom_divstab_agg_dat, iter = 5000, adapt_delta = 0.999)
dsra_gb_b <- stan_lmer(gamma_var ~ beta_div_mean + (beta_div_mean|lter_site), data = metacom_divstab_agg_dat, iter = 5000, adapt_delta = 0.99)
dsra_b_b  <- stan_lmer(phi_var ~ beta_div_mean + (beta_div_mean|lter_site), data = metacom_divstab_agg_dat, iter = 5000, adapt_delta = 0.99)
dsra_a_b  <- stan_lmer(alpha_var ~ alpha_div_mean + (alpha_div_mean|lter_site), data = metacom_divstab_agg_dat, iter = 5000, adapt_delta = 0.99)

# # nonlinear
# dsra_g_b_2  <- stan_lmer(gamma_var ~ gamma_div_mean + I(gamma_div_mean^2) + (gamma_div_mean|lter_site), data = metacom_divstab_agg_dat, iter = 5000, adapt_delta = 0.99)
# dsra_ga_b_2 <- stan_lmer(gamma_var ~ alpha_div_mean + I(alpha_div_mean^2) + (alpha_div_mean|lter_site), data = metacom_divstab_agg_dat, iter = 5000, adapt_delta = 0.99)
# dsra_gb_b_2 <- stan_lmer(gamma_var ~ beta_div_mean + I(beta_div_mean^2) + (beta_div_mean|lter_site), data = metacom_divstab_agg_dat, iter = 5000, adapt_delta = 0.99)
# dsra_b_b_2  <- stan_lmer(phi_var ~ beta_div_mean + I(beta_div_mean^2) + (beta_div_mean|lter_site), data = metacom_divstab_agg_dat, iter = 5000, adapt_delta = 0.99)
# dsra_a_b_2  <- stan_lmer(alpha_var ~ alpha_div_mean + I(alpha_div_mean^2) + (alpha_div_mean|lter_site), data = metacom_divstab_agg_dat, iter = 5000, adapt_delta = 0.99)

# # compare linear vs nonlinear models
# loo1 <- loo(dsrc_g_b, k_threshold = 0.7)
# loo2 <- loo(dsrc_g_b_2, k_threshold = 0.7)
# loo_compare(loo1, loo2) # quad better
# loo1 <- loo(dsrc_ga_b, k_threshold = 0.7)
# loo2 <- loo(dsrc_ga_b_2, k_threshold = 0.7)
# loo_compare(loo1, loo2) # linear better
# loo1 <- loo(dsrc_gb_b, k_threshold = 0.7)
# loo2 <- loo(dsrc_gb_b_2, k_threshold = 0.7)
# loo_compare(loo1, loo2) # linear better
# loo1 <- loo(dsrc_b_b, k_threshold = 0.7)
# loo2 <- loo(dsrc_b_b_2, k_threshold = 0.7)
# loo_compare(loo1, loo2) # quad better
# loo1 <- loo(dsrc_a_b, k_threshold = 0.7)
# loo2 <- loo(dsrc_a_b_2, k_threshold = 0.7)
# loo_compare(loo1, loo2)# linear better
# 
# loo1 <- loo(dsra_g_b, k_threshold = 0.7)
# loo2 <- loo(dsra_g_b_2, k_threshold = 0.7)
# loo_compare(loo1, loo2) # linear better
# loo1 <- loo(dsra_ga_b, k_threshold = 0.7)
# loo2 <- loo(dsra_ga_b_2, k_threshold = 0.7)
# loo_compare(loo1, loo2) # linear better
# loo1 <- loo(dsra_gb_b, k_threshold = 0.7)
# loo2 <- loo(dsra_gb_b_2, k_threshold = 0.7)
# loo_compare(loo1, loo2) # linear better
# loo1 <- loo(dsra_b_b, k_threshold = 0.7)
# loo2 <- loo(dsra_b_b_2, k_threshold = 0.7)
# loo_compare(loo1, loo2) # quad
# loo1 <- loo(dsra_a_b, k_threshold = 0.7)
# loo2 <- loo(dsra_a_b_2, k_threshold = 0.7)
# loo_compare(loo1, loo2) # linear 


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


# regional compositional variability

# compositional models

# gamma variability vs. gamma diversity
b0 <- tidy(dsrc_g_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsrc_g_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]

(div_stab_gamma_h <- metacom_divstab_comp_dat %>% 
    add_epred_draws(dsrc_g_b, ndraws = 200, re_formla = NA) %>% 
    ggplot(aes(x = gamma_div_mean, y = gamma_var, label = lter_site)) +
    geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
    geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
    labs(x = expression(paste("Mean ", gamma, "-diversity")),
         y = expression(paste("Comp. ", gamma, "-variability")), # (", BD[gamma]^h, ") for metric
         color = "Organism group", shape = "Organism group")  +
    scale_color_manual(values = pal, drop = FALSE)
)

# alpha variability vs. alpha diversity
b0 <- tidy(dsrc_a_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsrc_a_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]

(div_stab_alpha_h <- metacom_divstab_comp_dat %>% 
    add_epred_draws(dsrc_a_b, ndraws = 200, re_formla = NA) %>% 
    ggplot(aes(x = alpha_div_mean, y = alpha_var, label = lter_site)) +
    geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
    geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
    labs(x = expression(paste("Mean ", alpha, "-diversity")),
         y = expression(paste("Comp. ", alpha, "-variability")),
         color = "Organism group", shape = "Organism group") +
    scale_color_manual(values = pal, drop = FALSE)
)

# gamma variability vs. alpha diversity
b0 <- tidy(dsrc_ga_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsrc_ga_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]

(div_stab_alpha_gamma_h <- metacom_divstab_comp_dat %>% 
    add_epred_draws(dsrc_ga_b, ndraws = 200, re_formula = NA) %>% 
    ggplot(aes(x = alpha_div_mean, y = gamma_var, label = lter_site)) +
    geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
    geom_abline(intercept = b0, slope = b1, size = 1, color = "black") +
    geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
    labs(x = expression(paste("Mean ", alpha, "-diversity")),
         y = expression(paste("Comp. ", gamma, "-variability")),
         color = "Organism group", shape = "Organism group") +
    scale_color_manual(values = pal, drop = FALSE)
)

# gamma variability vs. beta diversity
b0 <- tidy(dsrc_gb_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsrc_gb_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]

(div_stab_beta_h <- metacom_divstab_comp_dat %>% 
    add_epred_draws(dsrc_gb_b, ndraws = 200, re_formula = NA) %>% 
    ggplot(aes(x = beta_div_mean, y = gamma_var, label = lter_site)) +
    geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
    geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
    scale_color_manual(values = pal, drop = FALSE) +
    labs(x = expression(paste("Mean ", beta, "-diversity")),
         y = expression(paste("Comp. ", gamma, "-variability")),
         color = "Organism group", shape = "Organism group")
)

# regional aggregate variability

# gamma variability vs. gamma diversity
b0 <- tidy(dsra_g_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsra_g_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]

(div_stab_gamma_agg <- metacom_divstab_agg_dat %>% 
    add_epred_draws(dsra_g_b, ndraws = 200, re_formula = NA) %>% 
    ggplot(aes(x = gamma_div_mean, y = gamma_var, label = lter_site)) +
    geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
    geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
    scale_color_manual(values = pal, drop = FALSE) +
    labs(x = expression(paste("Mean ", gamma, "-diversity")),
         y = expression(paste("Agg. ", gamma, "-variability")),
         color = "Organism group", shape = "Organism group") 
)

# alpha variability vs. alpha diversity
b0 <- tidy(dsra_a_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsra_a_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]

(div_stab_alpha_agg <- metacom_divstab_agg_dat %>% 
    add_epred_draws(dsra_a_b, ndraws = 200, re_formula = NA) %>% 
    ggplot(aes(x = alpha_div_mean, y = alpha_var, label = lter_site)) +
    geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
    geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
    scale_color_manual(values = pal, drop = FALSE) +
    labs(x = expression(paste("Mean ", alpha, "-diversity")),
         y = expression(paste("Agg. ", alpha, "-variability")),
         color = "Organism group", shape = "Organism group")
)

# gamma variability vs. alpha diversity
b0 <- tidy(dsra_ga_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsra_ga_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]

(div_stab_alpha_gamma_agg <- metacom_divstab_agg_dat %>% 
    add_epred_draws(dsra_ga_b, ndraws = 200, re_formula = NA) %>% 
    ggplot(aes(x = alpha_div_mean, y = gamma_var, label = lter_site)) +
    geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
    geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
    scale_color_manual(values = pal, drop = FALSE) +
    labs(x = expression(paste("Mean ", alpha, "-diversity")),
         y = expression(paste("Agg. ", gamma, "-variability")),
         color = "Organism group", shape = "Organism group")
)

# gamma variability vs. beta diversity
b0 <- tidy(dsra_gb_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsra_gb_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]

(div_stab_beta_agg <- metacom_divstab_agg_dat %>% 
    add_epred_draws(dsra_ga_b, ndraws = 200, re_formula = NA) %>% 
    ggplot(aes(x = beta_div_mean, y = gamma_var)) +
    geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
    geom_abline(intercept = b0, slope = b1, size = 1, color = "black") +
    geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
    scale_color_manual(values = pal, drop = FALSE) +
    labs(x = expression(paste("Mean ", beta, "-diversity")),
         y = expression(paste("Agg. ", gamma, "-variability")),
         color = "Organism group", shape = "Organism group")
)

# now combine panels into Figure 3
div_stab_multipanel_fig <- div_stab_gamma_agg + div_stab_gamma_h +
  div_stab_beta_agg + div_stab_beta_h +
  div_stab_alpha_gamma_agg + div_stab_alpha_gamma_h +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect", nrow = 3)
ggsave(filename = "figs/diversity_variability.png", plot = div_stab_multipanel_fig, width = 2.5*4, height = 3*3, units = "in", dpi = 1000, bg = "white")


# 4. COMPARE BETA DIVERSITY WITH PHI, SYNCHRONY

b0 <- tidy(dsrc_b_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsrc_b_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]

(div_stab_phi_h <- metacom_divstab_comp_dat %>% 
    add_epred_draws(dsrc_b_b, ndraws = 200, re_formula = NA) %>% 
    ggplot(aes(x = beta_div_mean, y = phi_var)) +
    geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
    geom_abline(intercept = b0, slope = b1, size = 1, color = "black") +
    geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
    scale_y_continuous(limits = c(0,1)) +
    scale_color_manual(values = pal, drop = FALSE) +
    labs(x = expression(paste("Mean ", beta, "-diversity")),
         y = expression(paste("Comp. Synchrony (", BD[phi]^h, ")")),
         color = "Organism group", shape = "Organism group")
)

b0 <- tidy(dsra_b_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[1]
b1 <- tidy(dsra_b_b, effects = "fixed", conf.int = TRUE, conf.level = 0.80)$estimate[2]

(div_stab_phi_agg <- metacom_divstab_agg_dat %>% 
    add_epred_draws(dsra_b_b, ndraws = 200, re_formula = NA) %>% 
    ggplot(aes(x = beta_div_mean, y = phi_var)) +
    geom_line(aes(y = .epred, group = .draw), color = "grey", alpha = 0.1) + 
    geom_abline(intercept = b0, slope = b1, size = 1, color = "black") +
    geom_point(size = 3, alpha = 1, mapping = aes(color = organism_group, shape = organism_group)) +
    scale_y_continuous(limits = c(0,1))+
    scale_color_manual(values = pal, drop = FALSE) +
    labs(x = expression(paste("Mean ", beta, "-diversity")),
         y = expression(paste("Agg. Synchrony (",phi,")")),
         color = "Organism group", shape = "Organism group") 
)


# combine panels to make Fig 1c and d
div_stab_phi_beta_fig <- (div_stab_phi_agg) + div_stab_phi_h + 
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", nrow = 2)
ggsave(filename = "figs/beta_phi.png", plot = div_stab_phi_beta_fig, width = 1.5*4, height = 2*3, units = "in", dpi = 1000, bg = "white")


# Combine all alpha-alpha, beta-beta panels to make Figure 1
div_stab_multiscale_partition_fig <- 
  (local_divstab_agg_fig + theme(legend.position = "none")) + (local_divstab_comp_fig + theme(legend.position = "none")) +
  div_stab_phi_agg + div_stab_phi_h +
  #div_stab_gamma_agg + div_stab_gamma_h +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect", nrow = 2)
ggsave(filename = "figs/diversity_variability_alpha_beta.png", plot = div_stab_multiscale_partition_fig, width = 2.7*4, height = 3.2*3, units = "in", dpi = 1000, bg = "white")



# 4. COMPOSITIONAL VS AGGREGATE VARIABILITY -----------------------------------
comp_agg_stab <- left_join(metacom_divstab_agg_dat, metacom_divstab_comp_dat, by = c("dataset_id", "lter_site", "organism_group"), suffix = c("_agg", "_comp"))

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
  annotate("text", x = .75, y = 0.1, size = 5, label = expression(paste(rho, "= 0.47")))
ggsave("figs/phi_comparison.png",plot = phi_compare, bg = "white", width = 6, height = 3/4*6, dpi = 600)


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


metacom_var_renamed <- str_split(metacom_var$dataset_id, pattern = "-", simplify = T)
colnames(metacom_var_renamed) <- c("site", "org", "person")
metacom_var_renamed <- metacom_var_renamed %>% as_tibble() %>% bind_cols(metacom_var) %>% 
  mutate(org = paste0(org, 
      ifelse(site == "sgs" & org == "plants" & person == "compagnoni", 2, 
      ifelse(site == "sgs" & org == "plants" & person == "catano", 1, "")))) %>% 
  mutate(dataset_id = paste0(site,"-",org))

values_by_dataset_fig <- metacom_var_renamed %>% 
  filter(metric %in% c("gamma_var", "alpha_var", "phi_var")) %>% 
  mutate(metric = ifelse(metric == "gamma_var", "Metacommunity",
                         ifelse(metric == "phi_var", "Spatial synchrony", "Local"))) %>% 
  mutate(metric = factor(metric, levels = c("Metacommunity", "Spatial synchrony", "Local"))) %>% 
  mutate(variability_type = ifelse(variability_type == "agg", 
                                   "Aggregate", "Compositional")) %>% 
  ggplot(aes(x = dataset_id, y = metric_value, label = round(metric_value, 2), color = organism_group)) + 
  geom_point() + 
  geom_text_repel(nudge_x = 0.1, show.legend = FALSE) +
  scale_color_manual(values = pal, drop = FALSE) +
  facet_grid(variability_type ~ metric) +
  coord_flip() +
  #theme(legend.position = "none") + 
  theme(panel.grid.major.x = element_line(color = "grey90")) +
  labs(x = "", y = "", color = "Organism group")
#values_by_dataset_fig
ggsave(plot = values_by_dataset_fig,
       filename = "figs/variability_dataset.png",
       width = 12, height = 12, dpi = 1000, bg = "white")
