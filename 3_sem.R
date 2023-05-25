library(piecewiseSEM)
library(DiagrammeR)
library(tidyverse)
library(lavaan)
library(lme4)
library(here)

env_var <- read_csv(here("data/lter_centroid_satdata.csv")) %>% 
  rename("lter_site" = site)
data_list <- read_csv(here("data/L3_DATA_list.csv")) %>% 
  filter(l3_filename != "L3-mcr-fish-castorani.csv") %>% 
  select(
    dataset_id, `LTER site`, mobility, `trophic group`, biome, organism_group, n.plots
  ) %>% 
  rename("lter_site" = "LTER site")
dsr_tot <- as.data.frame(read_csv(file = "dsr_results_table.csv"))

tot_mod <- psem(
  lm(cv_gamma ~ cv_alpha + phi + bd_alpha + bd_phi + bd_gamma + gamma_div_mean + alpha_div_mean + beta_div_mean, data = dsr_tot),
  lm(bd_gamma ~ bd_alpha + bd_phi + gamma_div_mean + alpha_div_mean + beta_div_mean, data = dsr_tot),
  lm(gamma_div_mean ~ alpha_div_mean + beta_div_mean, data = dsr_tot),
  lm(phi ~ beta_div_mean + bd_phi + cv_alpha + bd_alpha, data = dsr_tot),
  lm(cv_alpha ~ alpha_div_mean + bd_alpha + beta_div_mean + gamma_div_mean, data = dsr_tot),
  lm(bd_phi ~ beta_div_mean + bd_alpha, data = dsr_tot),
  lm(bd_alpha ~ alpha_div_mean + beta_div_mean + gamma_div_mean, data = dsr_tot),
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
         penwidth = abs(coefs(tot_mod)$Std.Estimate) * 5) 

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
sem_path_sigonly <- sem_path %>%
  select_edges(conditions = style != "solid") %>%
  delete_edges_ws()

sem_plot <- DiagrammeR::render_graph(sem_path)
sem_plot



sem_path_sigonly$edges_df <- sem_path_sigonly$edges_df %>% 
  mutate(style = ifelse(as.numeric(label) < 0, "dashed", "solid"))

sem_plot_sigonly <- DiagrammeR::render_graph(sem_path_sigonly)
sem_plot_sigonly



#############
dsr_env <- env_var %>%
  group_by(lter_site) %>%
  summarize(elev_sd_mean = mean(elevation_spatial_sd_slope_60km, na.rm = T),
            bathymetry_sd_mean = mean(bathymetry_spatial_sd_slope_60km, na.rm = T)) %>% 
  right_join(left_join(dsr_tot, data_list)) %>%
  filter(lter_site != "HAYS")

env_mod <- psem(
  lm(cv_gamma ~ cv_alpha + phi + bd_alpha + bd_phi + bd_gamma + gamma_div_mean + alpha_div_mean + beta_div_mean, data = dsr_env),
  lm(bd_gamma ~ bd_alpha + bd_phi + gamma_div_mean + alpha_div_mean + beta_div_mean, data = dsr_env),
  lm(gamma_div_mean ~ alpha_div_mean + beta_div_mean, data = dsr_env),
  lm(phi ~ beta_div_mean + bd_phi + cv_alpha + bd_alpha, data = dsr_env),
  lm(cv_alpha ~ alpha_div_mean + bd_alpha + beta_div_mean + gamma_div_mean, data = dsr_env),
  lm(bd_phi ~ beta_div_mean + bd_alpha, data = dsr_env),
  lm(bd_alpha ~ alpha_div_mean + beta_div_mean + gamma_div_mean, data = dsr_env),
  data = dsr_env
)
summary(env_mod)
AIC_psem(env_mod)

env_mod <- psem(
  lm(cv_gamma ~ cv_alpha + phi + bd_alpha + bd_phi + bd_gamma + gamma_div_mean + alpha_div_mean + beta_div_mean + n.plots + elev_sd_mean + bathymetry_sd_mean, data = dsr_env),
  lm(bd_gamma ~ bd_alpha + bd_phi + gamma_div_mean + alpha_div_mean + beta_div_mean + n.plots + elev_sd_mean + bathymetry_sd_mean, data = dsr_env),
  lm(gamma_div_mean ~ alpha_div_mean + beta_div_mean + n.plots + elev_sd_mean + bathymetry_sd_mean, data = dsr_env),
  lm(phi ~ beta_div_mean + bd_phi + cv_alpha + bd_alpha + n.plots + elev_sd_mean + bathymetry_sd_mean, data = dsr_env),
  lm(cv_alpha ~ alpha_div_mean + bd_alpha + beta_div_mean + gamma_div_mean + n.plots + elev_sd_mean + bathymetry_sd_mean, data = dsr_env),
  lm(bd_phi ~ beta_div_mean + bd_alpha + n.plots + elev_sd_mean + bathymetry_sd_mean, data = dsr_env),
  lm(bd_alpha ~ alpha_div_mean + beta_div_mean + gamma_div_mean + n.plots + elev_sd_mean + bathymetry_sd_mean, data = dsr_env),
  data = dsr_env
)
summary(env_mod)
AIC_psem(env_mod)

env_path <- plot(
  env_mod,
  return = T,
  #node_attrs = data.frame(shape = "rectangle", color = "black", fillcolor = "white"),
  #edge_attrs = data.frame(style = "solid", color = "black"),
  ns_dashed = T,
  alpha = 0.05,
  show = "std",
  digits = 3,
  add_edge_label_spaces = TRUE)

env_path$edges_df <- env_path$edges_df %>% 
  mutate(color = ifelse(as.numeric(label) < 0, "red", "blue"),
         penwidth = abs(coefs(env_mod)$Std.Estimate) * 5) 

env_path$edges_df

env_path$global_attrs[7, 2] <- "false"
env_path$global_attrs[1, 2] <- "dot"
env_path$nodes_df$label <- c(
  "Agg. Metacommunity \nVariability (CV_gamma)",
  "Comp. Metacommunity \nVariability (BD_gamma)",
  "Avg. Gamma Diversity",
  "Agg. Spatial \nSynchrony (phi)",
  "Agg. Local \nVariability (CV_alpha)",
  "Comp. Spatial \nSynchrony (BD_phi)",
  "Comp. Local \nVariability (BD_alpha)",
  "Avg. Alpha Diversity",
  "Avg. Beta Diversity",
  "Num. Plots", 
  "Elevation sd. mean",
  "Bathymetry sd. mean"
)
env_path$edges_df
env_path$global_attrs[5,2] <- 10
env_path$global_attrs[14,2] <- 12

# remove non significant edges
env_path_sigonly <- env_path %>%
  select_edges(conditions = style != "solid") %>%
  delete_edges_ws()

env_plot <- DiagrammeR::render_graph(env_path)
env_plot



env_path_sigonly$edges_df <- env_path_sigonly$edges_df %>% 
  mutate(style = ifelse(as.numeric(label) < 0, "dashed", "solid"))

env_plot_sigonly <- DiagrammeR::render_graph(env_path_sigonly)
env_plot_sigonly


