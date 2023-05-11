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
    dataset_id, `LTER site`, mobility, `trophic group`, biome, organism_group
  )
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

# indirect paths
0.991*0.089
0.594 * 0.089
-.62*.525

-.407*.683
.398*.835



# ##########
# mutate(dsr_tot, across(cv_gamma:bd_alpha, log10))
# comp_mod <- psem(
#   lm(bd_gamma ~ bd_phi + bd_alpha + gamma_div_mean + alpha_div_mean + beta_div_mean, data = dsr_tot),
#   lm(bd_phi ~ beta_div_mean, data = dsr_tot),
#   lm(bd_alpha ~ alpha_div_mean, data = dsr_tot),
#   data = dsr_tot
# )
# 
# summary(comp_mod)
# 
# 
# render_graph(plot(
#   comp_mod,
#   return = TRUE,
#   node_attrs = data.frame(shape = "rectangle", color = "black", fillcolor = "white"),
#   edge_attrs = data.frame(style = "solid", color = "black"),
#   ns_dashed = T,
#   alpha = 0.05,
#   show = "std",
#   digits = 3,
#   add_edge_label_spaces = TRUE), layout = "tree")
# 
# agg_mod <- psem(
#   lm(cv_gamma ~ phi + cv_alpha + gamma_div_mean, data = dsr_tot),
#   lm(phi ~ beta_div_mean, data = dsr_tot),
#   lm(cv_alpha ~ alpha_div_mean, data = dsr_tot),
#   data = dsr_tot
# )
# 
# summary(agg_mod)
# render_graph(plot(agg_mod, return = TRUE), layout = "tree")
# 
# 
# 
# #####
# 
# dsr_tot <- dsr_tot %>% left_join(data_list) %>% 
#   rename(lter_site = `LTER site`)
# 
# hist(log(dsr_tot$cv_alpha))
# hist(sqrt(dsr_tot$cv_gamma))
# hist(dsr_tot$phi)
# hist(log(dsr_tot$bd_alpha))
# hist(log(dsr_tot$bd_gamma))
# hist(dsr_tot$bd_phi)
# hist(log(dsr_tot$gamma_div_mean))
# hist(log(dsr_tot$beta_div_mean))
# hist(log(dsr_tot$alpha_div_mean))
# 
# dsr_tot_transformed <- dsr_tot %>% 
#   mutate(
#     cv_alpha = log(cv_alpha),
#     cv_gamma = sqrt(cv_gamma),
#     bd_alpha = log(bd_alpha),
#     bd_gamma = log(bd_gamma),
#     gamma_div_mean = log(gamma_div_mean),
#     beta_div_mean = log(beta_div_mean),
#     alpha_div_mean = log(alpha_div_mean)
#   )
# 
# summary(glm(cv_gamma ~ cv_alpha + phi + bd_alpha + bd_phi + bd_gamma + gamma_div_mean, family = Gamma, data = dsr_tot))
# summary(lm(cv_gamma ~ cv_alpha + phi + bd_alpha + bd_phi + bd_gamma + gamma_div_mean, data = dsr_tot))
# 
# g1 <- glm(cv_gamma ~ cv_alpha + phi + bd_alpha + bd_phi + bd_gamma + gamma_div_mean, family = gaussian(link = "log"), data = dsr_tot)
# step(g1, test = "LRT")
# glm(formula = cv_gamma ~ cv_alpha + phi + bd_alpha, family = gaussian(link = "log"), 
#     data = dsr_tot)
# 
# g2 <- glm(bd_gamma ~ bd_alpha + bd_phi + gamma_div_mean, family = gaussian(link = "log"), data = dsr_tot)
# step(g2, test = "LRT")
# 
# tot_mod <- psem(
#   lm(cv_gamma ~ cv_alpha + phi + bd_gamma, data = dsr_tot),
#   lm(bd_gamma ~ bd_alpha + bd_phi + gamma_div_mean, data = dsr_tot),
#   lm(gamma_div_mean ~ alpha_div_mean + beta_div_mean, data = dsr_tot),
#   lm(phi ~ beta_div_mean + bd_phi, data = dsr_tot),
#   lm(cv_alpha ~ alpha_div_mean + bd_alpha, data = dsr_tot),
#   lm(bd_phi ~ beta_div_mean, data = dsr_tot),
#   lm(bd_alpha ~ alpha_div_mean, data = dsr_tot),
#   data = dsr_tot
# )
# summary(tot_mod)
# 
# 
# tot_mod <- psem(
#   lm(cv_gamma ~ cv_alpha + phi, data = dsr_tot),
#   lm(bd_gamma ~ bd_alpha + bd_phi + gamma_div_mean, data = dsr_tot),
#   lm(gamma_div_mean ~ alpha_div_mean + beta_div_mean, data = dsr_tot),
#   lm(phi ~ beta_div_mean, data = dsr_tot),
#   lm(cv_alpha ~ bd_alpha, data = dsr_tot),
#   lm(bd_phi ~ beta_div_mean, data = dsr_tot),
#   data = dsr_tot
# )
# summary(tot_mod)
# 
# sem_path <- plot(
#   tot_mod,
#   return = T,
#   #node_attrs = data.frame(shape = "rectangle", color = "black", fillcolor = "white"),
#   #edge_attrs = data.frame(style = "solid", color = "black"),
#   ns_dashed = T,
#   alpha = 0.05,
#   show = "std",
#   digits = 3,
#   add_edge_label_spaces = TRUE)
# 
# sem_path$edges_df <- sem_path$edges_df %>% 
#   mutate(color = ifelse(as.numeric(label) < 0, "red", "blue"),
#          penwidth = abs(coefs(tot_mod)$Std.Estimate) * 10)
# 
# sem_path$edges_df
# 
# get_global_graph_attr_info(sem_path)
# sem_path$global_attrs[7, 2] <- "false"
# sem_path$global_attrs[1, 2] <- "dot"
# ?layout_nodes_w_string
# sem_path$nodes_df$label <- c(
#   "Agg. Metacommunity \nVariability (CV_gamma)",
#   "Comp. Metacommunity \nVariability (BD_gamma)",
#   "Avg. Gamma Diversity",
#   "Agg. Spatial \nSynchrony (phi)",
#   "Agg. Local \nVariability (CV_alpha)",
#   "Comp. Spatial \nSynchrony (BD_phi)",
#   "Comp. Local \nVariability (BD_alpha)",
#   "Avg. Alpha Diversity",
#   "Avg. Beta Diversity"
# )
# sem_path$edges_df
# sem_path$global_attrs[5,2] <- 10
# sem_path$global_attrs[14,2] <- 12
# 
# # remove non significant edges
# sem_path_sigonly <- sem_path %>%
#   select_edges(conditions = style != "solid") %>%
#   delete_edges_ws()
# 
# sem_plot <- DiagrammeR::render_graph(sem_path)
# sem_plot
# 
# 
# sem_plot <- DiagrammeR::render_graph(sem_path_sigonly)
# sem_plot
# 
# # indirect paths
# 0.991*0.089
# 0.594 * 0.089
# -.62*.525
# 
# -.407*.683
# .398*.835
# 
# 
# 
# #############
# dsr_tot <- env_var %>% 
#   group_by(lter_site) %>% 
#   summarize(elev_sd_mean = mean(elevation_spatial_sd_slope_60km),
#             bathymetry_sd_mean = mean(bathymetry_spatial_sd_slope_60km)) %>% print(n=100) %>% 
#   right_join(dsr_tot) %>% 
#   filter(lter_site != "HAYS")
# 
# tot_mod <- psem(
#   lm(cv_gamma ~ cv_alpha + phi, data = dsr_tot),
#   lm(bd_gamma ~ bd_alpha + bd_phi + gamma_div_mean, data = dsr_tot),
#   lm(gamma_div_mean ~ alpha_div_mean + beta_div_mean, data = dsr_tot),
#   lm(phi ~ beta_div_mean + elev_sd_mean + bathymetry_sd_mean, data = dsr_tot),
#   lm(cv_alpha ~ bd_alpha, data = dsr_tot),
#   lm(bd_phi ~ beta_div_mean + elev_sd_mean + bathymetry_sd_mean, data = dsr_tot),
#   lm(beta_div_mean ~ elev_sd_mean + bathymetry_sd_mean, data = dsr_tot),
#   data = dsr_tot
# )
# summary(tot_mod)


## 
# beta-div on agg mc var
-.718*.456*.651 + 1.124*(.407*.456*.651 + .49*.829) + .594*-1.274*(.407*.456*.651+.49*.829)

# beta-div on comp mc var
-.718*.496 + 1.124*(.407*.496+.721) + .594*-1.274*(.407*.496+.721)
