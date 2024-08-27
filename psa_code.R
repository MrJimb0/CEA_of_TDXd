library(dampack)  # To use CEA and PSA visualization 
library(patchwork)  # For combining ggplot2 figures

df_psa_lng <- read.csv("data/base_case_output.csv")
# Transform long to wide the dataset with group as columns
df_psa <- reshape(df_psa_lng[, -1], idvar = "sim", timevar = "group", direction = "wide")
# v_names_str <- df_psa_lng[1:4, "group"]
v_names_str <- c("Chemo → Chemo",
                 "Chemo → TDXd", 
                 "TDXd → Chemo", 
                 "TDXd → SG")
# expression(paste("Chemo", ""%right%"", "Chemo", sep = ""))
df_costs <- df_psa[, c(2, 5, 8, 11)]
df_qalys <- df_psa[, c(3, 6, 9, 12)]
df_lys   <- df_psa[, c(4, 7, 10, 13)]

## Visualize PSA results for CEA ----
### Create PSA object ----
l_psa <- dampack::make_psa_obj(cost          = df_costs, 
                               effectiveness = df_qalys, 
                               strategies    = v_names_str)
l_psa$strategies <- v_names_str
colnames(l_psa$effectiveness) <- v_names_str
colnames(l_psa$cost) <- v_names_str

#* Vector with willingness-to-pay (WTP) thresholds.
v_wtp <- seq(0, 700000, by = 5000)

### Cost-Effectiveness Scatter plot ----
txtsize <- 16
gg_scattter <- plot(l_psa, txtsize = txtsize) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  scale_y_continuous("Cost (Thousand $)", 
                     n.breaks = 6,
                     labels = function(x) x/1000) +
  xlab("Effectiveness (QALYs)") +
  guides(col = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom")
gg_scattter
ggsave(gg_scattter, filename = "figs/cea_scatter.png", width = 8, height = 6, dpi = 300)

### Incremental cost-effectiveness ratios (ICERs) with probabilistic output ----
#* Compute expected costs and effects for each strategy from the PSA
df_out_ce_psa <- summary(l_psa)
df_cea_psa <- dampack::calculate_icers(cost       = df_out_ce_psa$meanCost, 
                                       effect     = df_out_ce_psa$meanEffect,
                                       strategies = df_out_ce_psa$Strategy)
df_cea_psa

### Plot cost-effectiveness frontier with probabilistic output ----
plot(df_cea_psa, label = "all", txtsize = txtsize) +
  theme(legend.position = c(0.8, 0.2))

### Cost-effectiveness acceptability curves (CEACs) and frontier (CEAF) ---
ceac_obj <- dampack::ceac(wtp = v_wtp, psa = l_psa)
#* Regions of highest probability of cost-effectiveness for each strategy
summary(ceac_obj)
#* CEAC & CEAF plot
gg_ceac <- plot(ceac_obj, 
                txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 16) +
  xlab("Cost-effectiveness threshold (Thousand $/QALY)") +
  ylab("Probability of being cost-effective") +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme(legend.position = c(0.2, 0.48))
gg_ceac
ggsave(gg_ceac, filename = "figs/ceac.png", width = 8, height = 6, dpi = 300)

### Expected Loss Curves (ELCs) ----
elc_obj <- dampack::calc_exp_loss(wtp = v_wtp, psa = l_psa)
elc_obj

#* ELC plot
gg_elc <- plot(elc_obj, log_y = FALSE, 
               txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14,
               col = "full") +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  # geom_point(aes(shape = as.name("Strategy"))) +
  scale_y_continuous("Expected Loss (Thousand $)", 
                     n.breaks = 10,
                     labels = function(x) x/1000) +
  theme(legend.position = c(0.4, 0.7),)
gg_elc

### Expected value of perfect information (EVPI) ----
#* Function included in "R/Functions.R". The latest version can be found in `dampack` package
evpi <- dampack::calc_evpi(wtp = v_wtp, psa = l_psa)
#* EVPI plot
gg_evpi <- plot(evpi, effect_units = "QALY", 
                     txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14) +
  scale_y_continuous("EVPI (Thousand $)", 
                     n.breaks = 10,
                     labels = function(x) x/1000)
gg_evpi

### Combine all figures into one ----
patched_cea <- (gg_scattter +  gg_ceac + plot_layout(guides = "keep"))/(gg_elc + gg_evpi)
gg_psa_plots <- patched_cea + 
  plot_annotation(tag_levels = 'A')
gg_psa_plots

# Pairwise comparisons ----
## Chemo-Chemo vs Chemo-TDxD ----
l_psa_chemochemo_chemotdxd <- dampack::make_psa_obj(cost          = df_costs[, c(1, 2)],
                                                    effectiveness = df_qalys[, c(1, 2)],
                                                    strategies    = v_names_str[c(1, 2)])
### Incremental cost-effectiveness ratios (ICERs) with probabilistic output ----
#* Compute expected costs and effects for each strategy from the PSA
df_out_ce_psa_chemochemo_chemotdxd <- summary(l_psa_chemochemo_chemotdxd)
df_cea_psa_chemochemo_chemotdxd <- dampack::calculate_icers(cost       = df_out_ce_psa_chemochemo_chemotdxd$meanCost, 
                                       effect     = df_out_ce_psa_chemochemo_chemotdxd$meanEffect,
                                       strategies = df_out_ce_psa_chemochemo_chemotdxd$Strategy)
df_cea_psa_chemochemo_chemotdxd
### Plot cost-effectiveness frontier with probabilistic output ----
plot(df_cea_psa_chemochemo_chemotdxd, label = "all", txtsize = txtsize) +
  theme(legend.position = c(0.8, 0.2))

### Cost-effectiveness acceptability curves (CEACs) and frontier (CEAF) ---
ceac_obj_chemochemo_chemotdxd <- dampack::ceac(wtp = v_wtp, psa = l_psa_chemochemo_chemotdxd)
#* Regions of highest probability of cost-effectiveness for each strategy
summary(ceac_obj_chemochemo_chemotdxd)
#* CEAC & CEAF plot
gg_ceac_chemochemo_chemotdxd <- plot(ceac_obj_chemochemo_chemotdxd, txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme(legend.position = c(0.8, 0.48))
gg_ceac_chemochemo_chemotdxd

## Chemo-Chemo vs TDxd-Chemo ----
l_psa_chemotdxd_chemochemo <- dampack::make_psa_obj(cost          = df_costs[, c(1, 3)],
                                                    effectiveness = df_qalys[, c(1, 3)],
                                                    strategies    = v_names_str[c(1, 3)])
