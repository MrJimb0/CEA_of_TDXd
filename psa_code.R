library(dampack)  # To use CEA and PSA visualization 
library(patchwork)  # For combining ggplot2 figures

df_psa_lng <- read.csv("data/base_case_output.csv")
# Transform long to wide the dataset with group as columns
df_psa <- reshape(df_psa_lng[, -1], idvar = "sim", timevar = "group", direction = "wide")
v_names_str <- df_psa_lng[1:4, "group"]
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
txtsize <- 14
gg_scattter <- plot(l_psa, txtsize = txtsize) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  scale_y_continuous("Cost (Thousand $)",
                     labels = function(x) x/1000) +
  xlab("Effectiveness (QALYs)") +
  guides(col = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom")
gg_scattter

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
gg_ceac <- plot(ceac_obj, txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme(legend.position = c(0.8, 0.48))
gg_ceac

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
