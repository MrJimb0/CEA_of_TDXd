#This code was created by Marcus Moen, MS and edited by James Dickerson, MD, MS and Wesley Suen
#Oversight from Jeremy Goldhaber-Feibert, PhD and Fernando Alarid-Escudero, Ph.D.

options(scipen=999)
setwd("/Users/jamesdickerson/Library/CloudStorage/Box-Box/Dickerson Lab/Dickerson_Lab_Github/CEA_of_TDXd/")

#Plot file
library(ggplot2)
library(ggpubr)
library(dampack)
library(grid)
library(gridExtra)

source("finding_transition_probabilities_7state.R")
source("finding_transition_probabilities_8state.R")

#Store the transition matrices from the "finding_transition_probabilities_8state.R" file
A_tdxd_chemo <- tm_tdxd_chemo
A_chemo_chemo <- tm_chemo_chemo
A_chemo_tdxd <- tm_chemo_tdxd
A_tdxd_sg <- tm_tdxd_sg

plot_function <- function(input_var){
  A_tdxd_chemo <- input_var[[1]]
  A_chemo_chemo <- input_var[[2]]
  A_chemo_tdxd <- input_var[[3]]
  A_tdxd_sg <- input_var[[4]]
  
  
  #The extracted Kaplan-Meier values
  km_os_chemo_chemo <- c(1, 0.986, 0.981, 0.957, 0.94, 0.927, 0.884, 0.841, 0.793, 0.738, 0.712, 0.683, 0.669, 0.634)#, 0.610, 0.566, 0.520, 0.484, 0.459, 0.460, 0.457, 0.419, 0.375)
  km_pf_chemo_chemo <- c(1, 0.983, 0.756, 0.62, 0.603, 0.501, 0.439, 0.398, 0.34, 0.277, 0.268, 0.226, 0.217, 0.185) 
  km_os_tdxd_chemo <- c(1, 0.992, 0.985, 0.967, 0.958, 0.939, 0.928, 0.899, 0.87, 0.857, 0.829, 0.815, 0.793, 0.774)#, 0.742, 0.731, 0.688, 0.661, 0.623, 0.586, 0.566, 0.539, 0.507, 0.500, 0.484)
  km_pf_tdxd_chemo <- c(1, 0.988, 0.896, 0.82, 0.809, 0.765, 0.681, 0.633, 0.592, 0.554, 0.491, 0.457, 0.422, 0.386)#, 0.371, 0.363, 0.339, 0.31, 0.295, 0.274, 0.259) 
  
  #Vectors for storing the estimates of the KM curves
  km_os_chemo_chemo_model <- c()
  km_pf_chemo_chemo_model <- c()
  
  km_os_tdxd_chemo_model <- c()
  km_pf_tdxd_chemo_model <- c() 
  
  km_os_chemo_tdxd_model <- c()
  km_pf_chemo_tdxd_model <- c() 
  
  km_os_tdxd_sg_model <- c()
  km_pf_tdxd_sg_model <- c() 
  
  
  n <- 24
  #n = max(length(km_os_chemo_chemo), length(km_pf_chemo_chemo), length(km_os_tdxd_chemo), length(km_pf_tdxd_chemo))
  
  #The initial state
  s0_chemo_chemo <- matrix(c(1,0,0,0,0,0,0,0), ncol = 8)
  s0_tdxd_chemo <- matrix(c(1,0,0,0,0,0,0,0), ncol = 8)
  s0_chemo_tdxd <- matrix(c(1,0,0,0,0,0,0,0), ncol = 8)
  s0_tdxd_sg <- matrix(c(1,0,0,0,0,0,0,0), ncol = 8)
  
  survival_chemo_chemo <- s0_chemo_chemo[1]+s0_chemo_chemo[4]+s0_chemo_chemo[5]
  km_os_chemo_chemo_model <- append(km_os_chemo_chemo_model, survival_chemo_chemo)
  km_pf_chemo_chemo_model <- append(km_pf_chemo_chemo_model, s0_chemo_chemo[1])
  
  survival_tdxd_chemo <- s0_tdxd_chemo[1]+s0_tdxd_chemo[4]+s0_tdxd_chemo[5]
  km_os_tdxd_chemo_model <- append(km_os_tdxd_chemo_model, survival_tdxd_chemo)
  km_pf_tdxd_chemo_model <- append(km_pf_tdxd_chemo_model, s0_tdxd_chemo[1])
  
  survival_chemo_tdxd <- s0_chemo_tdxd[1]+s0_chemo_tdxd[4]+s0_chemo_tdxd[5]
  km_os_chemo_tdxd_model <- append(km_os_chemo_tdxd_model, survival_chemo_tdxd)
  km_pf_chemo_tdxd_model <- append(km_pf_chemo_tdxd_model, s0_chemo_tdxd[1])
  
  survival_tdxd_sg <- s0_tdxd_sg[1]+s0_tdxd_sg[4]+s0_tdxd_sg[5]
  km_os_tdxd_sg_model <- append(km_os_tdxd_sg_model, survival_tdxd_sg)
  km_pf_tdxd_sg_model <- append(km_pf_tdxd_sg_model, s0_tdxd_sg[1])
  
  #ae_test_chemo <- c(0)
  #ae_test_tdxd <- c(0)
  #ae_test_ild <- c(0)
  
  
  for(t in 1:n){
    
    s1_chemo_chemo <- s0_chemo_chemo %*% A_chemo_chemo[[t]]
    s1_tdxd_chemo <- s0_tdxd_chemo %*% A_tdxd_chemo[[t]]
    s1_chemo_tdxd <- s0_chemo_tdxd %*% A_chemo_tdxd[[t]]
    s1_tdxd_sg <- s0_tdxd_sg %*% A_tdxd_sg[[t]]
    
    survival_chemo_chemo <- km_os_chemo_chemo_model[t] * (s1_chemo_chemo[1] + s1_chemo_chemo[4] + s1_chemo_chemo[5]) / 
      (s1_chemo_chemo[1] + s1_chemo_chemo[4] + s1_chemo_chemo[5] + s0_chemo_chemo[1] * A_chemo_chemo[[t]][1, 8] + 
         s0_chemo_chemo[4] * A_chemo_chemo[[t]][4, 8] + s0_chemo_chemo[5] * A_chemo_chemo[[t]][5, 8])
    km_os_chemo_chemo_model <- append(km_os_chemo_chemo_model, survival_chemo_chemo)
    km_pf_chemo_chemo_model <- append(km_pf_chemo_chemo_model, km_pf_chemo_chemo_model[t] * s1_chemo_chemo[1] / 
                                        (s1_chemo_chemo[1] + s0_chemo_chemo[1] * A_chemo_chemo[[t]][1, 4] + 
                                           s0_chemo_chemo[1] * A_chemo_chemo[[t]][1, 5] + s0_chemo_chemo[1] * A_chemo_chemo[[t]][1, 8]))
    
    survival_tdxd_chemo <- km_os_tdxd_chemo_model[t] * (s1_tdxd_chemo[1] + s1_tdxd_chemo[4] + s1_tdxd_chemo[5]) / 
      (s1_tdxd_chemo[1] + s1_tdxd_chemo[4] + s1_tdxd_chemo[5] + s0_tdxd_chemo[1] * A_tdxd_chemo[[t]][1, 8] + 
         s0_tdxd_chemo[4] * A_tdxd_chemo[[t]][4, 8] + s0_tdxd_chemo[5] * A_tdxd_chemo[[t]][5, 8])
    km_os_tdxd_chemo_model <- append(km_os_tdxd_chemo_model, survival_tdxd_chemo)
    km_pf_tdxd_chemo_model <- append(km_pf_tdxd_chemo_model, km_pf_tdxd_chemo_model[t] * s1_tdxd_chemo[1] / 
                                       (s1_tdxd_chemo[1] + s0_tdxd_chemo[1] * A_tdxd_chemo[[t]][1, 4] + 
                                          s0_tdxd_chemo[1] * A_tdxd_chemo[[t]][1, 5] + s0_tdxd_chemo[1] * A_tdxd_chemo[[t]][1, 8]))
    
    survival_chemo_tdxd <- km_os_chemo_tdxd_model[t] * (s1_chemo_tdxd[1] + s1_chemo_tdxd[4] + s1_chemo_tdxd[5]) / 
      (s1_chemo_tdxd[1] + s1_chemo_tdxd[4] + s1_chemo_tdxd[5] + s0_chemo_tdxd[1] * A_chemo_tdxd[[t]][1, 8] + 
         s0_chemo_tdxd[4] * A_chemo_tdxd[[t]][4, 8] + s0_chemo_tdxd[5] * A_chemo_tdxd[[t]][5, 8])
    km_os_chemo_tdxd_model <- append(km_os_chemo_tdxd_model, survival_chemo_tdxd)
    km_pf_chemo_tdxd_model <- append(km_pf_chemo_tdxd_model, km_pf_chemo_tdxd_model[t] * s1_chemo_tdxd[1] / 
                                       (s1_chemo_tdxd[1] + s0_chemo_tdxd[1] * A_chemo_tdxd[[t]][1, 4] + 
                                          s0_chemo_tdxd[1] * A_chemo_tdxd[[t]][1, 5] + s0_chemo_tdxd[1] * A_chemo_tdxd[[t]][1, 8]))
    
    survival_tdxd_sg <- km_os_tdxd_sg_model[t] * (s1_tdxd_sg[1] + s1_tdxd_sg[4] + s1_tdxd_sg[5]) / 
      (s1_tdxd_sg[1] + s1_tdxd_sg[4] + s1_tdxd_sg[5] + s0_tdxd_sg[1] * A_tdxd_sg[[t]][1, 8] + 
         s0_tdxd_sg[4] * A_tdxd_sg[[t]][4, 8] + s0_tdxd_sg[5] * A_tdxd_sg[[t]][5, 8])
    km_os_tdxd_sg_model <- append(km_os_tdxd_sg_model, survival_tdxd_sg)
    km_pf_tdxd_sg_model <- append(km_pf_tdxd_sg_model, km_pf_tdxd_sg_model[t] * s1_tdxd_sg[1] / 
                                    (s1_tdxd_sg[1] + s0_tdxd_sg[1] * A_tdxd_sg[[t]][1, 4] + 
                                       s0_tdxd_sg[1] * A_tdxd_sg[[t]][1, 5] + s0_tdxd_sg[1] * A_tdxd_sg[[t]][1, 8]))
    
    s0_chemo_chemo <- s1_chemo_chemo
    s0_tdxd_chemo <- s1_tdxd_chemo
    s0_chemo_tdxd <- s1_chemo_tdxd
    s0_tdxd_sg <- s1_tdxd_sg
    
  }
  
  m = 24
  #m = min(length(km_os_chemo_chemo), length(km_pf_chemo_chemo), length(km_os_tdxd_chemo), length(km_pf_tdxd_chemo))
  idx <- 1:m
  
  # Create the data frame for plotting
  df <- data.frame(
    idx = idx,
    km_os_chemo_chemo_model = km_os_chemo_chemo_model[1:m],
    km_os_tdxd_chemo_model = km_os_tdxd_chemo_model[1:m],
    km_os_chemo_chemo = km_os_chemo_chemo[1:m],
    km_os_tdxd_chemo = km_os_tdxd_chemo[1:m],
    km_os_chemo_tdxd_model = km_os_chemo_tdxd_model[1:m],
    km_os_tdxd_sg_model = km_os_tdxd_sg_model[1:m]
  )
  
  # Convert the y-values to percentages
  df$km_os_chemo_chemo_model <- df$km_os_chemo_chemo_model * 100
  df$km_os_tdxd_chemo_model <- df$km_os_tdxd_chemo_model * 100
  df$km_os_chemo_chemo <- df$km_os_chemo_chemo * 100
  df$km_os_tdxd_chemo <- df$km_os_tdxd_chemo * 100
  df$km_os_chemo_tdxd_model <- df$km_os_chemo_tdxd_model * 100
  df$km_os_tdxd_sg_model <- df$km_os_tdxd_sg_model * 100  
  
  plot1 <- ggplot(df, aes(x = idx)) + 
    geom_line(aes(y = km_os_tdxd_chemo, color = "T-DXd->Chemo (KM)"), linetype = "twodash", size = 1.2) +
    geom_line(aes(y = km_os_tdxd_sg_model, color = "T-DXd->SG (modeled)"), size = 1.2) +
    geom_line(aes(y = km_os_tdxd_chemo_model, color = "T-DXd->Chemo (modeled)"), size = 1.2) +
    geom_line(aes(y = km_os_chemo_tdxd_model, color = "Chemo->T-DXd (modeled)"), size = 1.2) +
    geom_line(aes(y = km_os_chemo_chemo_model, color = "Chemo Only (modeled)"), size = 1.2) +
    geom_line(aes(y = km_os_chemo_chemo, color = "Chemo Only (KM)"), linetype = "twodash", size = 1.2) +
    scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100), expand = c(0, 0)) +
    ylab("OS Probability (%)") +
    xlab("Months") +
    ggtitle("Modeled OS vs. Destiny-Breast04 Trial KM Estimates") +
    scale_x_continuous(breaks = round(seq(0, 24, by = 3), 1), limits = c(0, 24), expand = c(0, 0)) +
    theme_classic(base_size = 12) + 
    theme(
      text = element_text(family = "Arial"),
      axis.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "right",  # Place legend at the bottom of the plot
      legend.text = element_text(family = "Arial", size = 12),  # Set legend text font and size
      legend.title = element_text(family = "Arial", size = 14, face = "bold")  # Set legend title font and size
    ) +
    scale_color_manual(
      values = c(
        "Chemo Only (modeled)" = "red",
        "Chemo Only (KM)" = "blue",
        "T-DXd->Chemo (modeled)" = "chocolate",
        "T-DXd->Chemo (KM)" = "green",
        "Chemo->T-DXd (modeled)" = "deeppink",
        "T-DXd->SG (modeled)" = "darkkhaki"
      ),
      breaks = c(
        "T-DXd->Chemo (KM)",
        "T-DXd->SG (modeled)",
        "T-DXd->Chemo (modeled)",
        "Chemo->T-DXd (modeled)",
        "Chemo Only (modeled)",
        "Chemo Only (KM)"
      )
    )
  
  # Add manual legend
  plot1 <- plot1 + guides(color = guide_legend(title = "Treatment Strategies"))  # Add legend title
  
  print(plot1)
  
  
  
  plot2 <- ggplot(df, aes(x = idx)) + 
    geom_line(aes(y = km_pf_chemo_chemo_model[1:m]), color = "red") + 
    geom_line(aes(y = km_pf_chemo_chemo[1:m]), color = "blue", linetype = "twodash") +
    geom_line(aes(y = km_pf_tdxd_chemo_model[1:m]), color = "chocolate") + 
    geom_line(aes(y = km_pf_tdxd_chemo[1:m]), color = "green", linetype = "twodash") +
    geom_line(aes(y = km_pf_chemo_tdxd_model[1:m]), color = "deeppink") + 
    geom_line(aes(y = km_pf_tdxd_sg_model[1:m]), color = "darkkhaki") + 
    ylim(0, 1) +
    ylab("Probability") +
    xlab("Month") +
    ggtitle("PF") +
    scale_x_continuous(breaks = round(seq(0, 18, by = 3), 1)) +
    theme_minimal() + 
    theme(
      text = element_text(family = "Arial"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 16)
    )
  
  print(plot2)
}

opt_var <- list(A_tdxd_chemo,
                A_chemo_chemo,
                A_chemo_tdxd,
                A_tdxd_sg)

plot_function(opt_var)














source("finding_transition_probabilities_8state.R")

plot_function <- function(tm_tdxd_chemo){
  
  
  A_list <- tm_tdxd_chemo
  
  
  #The extracted Kaplan-Meier values
  km_os_tdxd <- c(1, 0.992, 0.985, 0.967, 0.958, 0.939, 0.928, 0.899, 0.87, 0.857, 0.829, 0.815, 0.793, 0.774)#, 0.742, 0.731, 0.688, 0.661, 0.623, 0.586, 0.566, 0.539, 0.507, 0.500, 0.484)
  km_pf_tdxd <- c(1, 0.988, 0.896, 0.82, 0.809, 0.765, 0.681, 0.633, 0.592, 0.554, 0.491, 0.457, 0.422, 0.386)#, 0.371, 0.363, 0.339, 0.31, 0.295, 0.274, 0.259) 
  
  #Vectors for storing the estimates of the KM curves
  km_os_tdxd_model <- c()
  km_pf_tdxd_model <- c() 
  
  
  n = max(length(km_os_tdxd), length(km_pf_tdxd))
  
  #Start states
  s0_tdxd <- c(1,0,0,0,0,0,0,0)
  
  survival_tdxd <- s0_tdxd[1]+s0_tdxd[4]+s0_tdxd[5]
  km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
  km_pf_tdxd_model <- append(km_pf_tdxd_model, s0_tdxd[1])
  
  ae_test_tdxd <- c(0)
  ae_test_ild <- c(0)
  
  
  for(t in 1:n){
    A_tdxd <- A_list[[t]]
    s1_tdxd <- s0_tdxd %*% A_tdxd
    
    
    survival_tdxd <- km_os_tdxd_model[t]*(s1_tdxd[1]+s1_tdxd[4]+s1_tdxd[5])/(s1_tdxd[1]+s1_tdxd[4]+s1_tdxd[5]+s0_tdxd[1]*A_tdxd[1,8]+s0_tdxd[4]*A_tdxd[4,8]+s0_tdxd[5]*A_tdxd[5,8])
    km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
    km_pf_tdxd_model <- append(km_pf_tdxd_model, km_pf_tdxd_model[t]*(s1_tdxd[1])/(s1_tdxd[1]+s0_tdxd[1]*A_tdxd[1,4]+s0_tdxd[1]*A_tdxd[1,5]+s0_tdxd[1]*A_tdxd[1,8]))
    
    ae_test_tdxd <- c(ae_test_tdxd, ae_test_tdxd[t]+s0_tdxd[1]*A_tdxd[1,2])
    ae_test_ild <- c(ae_test_ild, ae_test_ild[t]+s0_tdxd[1]*A_tdxd[1,3])
    
    s0_tdxd <- s1_tdxd
    
  }
  
  
  
  m = min(length(km_os_tdxd), length(km_pf_tdxd))
  idx <- 1:m
  df <- data.frame(idx, km_os_tdxd_model[1:m])
  
  
  
  plot2 <- ggplot(df, aes(x=idx)) + 
    geom_line(aes(y = km_os_tdxd_model[1:m]), color = "red") + 
    geom_line(aes(y = km_os_tdxd[1:m]), color="blue", linetype="twodash") +
    ylim(0, 1) +
    ylab("Probability")+
    xlab("Month")+
    ggtitle("OS T-DxD")+
    scale_x_continuous(breaks = round(seq(0, 18, by = 3),1)) + 
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))
  
  plot4 <- ggplot(df, aes(x = idx)) + 
    geom_line(aes(y = km_pf_tdxd_model[1:m]), color = "red") + 
    geom_line(aes(y = km_pf_tdxd[1:m]), color = "blue", linetype = "twodash") +
    ylim(0, 1) +
    ylab("Probability") +
    xlab("Month") +
    ggtitle("PFS T-DxD") +
    scale_x_continuous(breaks = round(seq(0, 18, by = 3),1)) + 
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))
  
  print(ae_test_tdxd)
  print(ae_test_ild)
  plot(0:(length(ae_test_tdxd)-1), ae_test_tdxd, type = "l")
  plot(0:(length(ae_test_tdxd)-1), ae_test_ild, type = "l")
  
  ggarrange(plot2, plot4,
            ncol = 1, nrow = 2, common.legend = TRUE,legend="bottom") 
  
}

plot_function(tm_tdxd_chemo)


















source("CEA_8state.R")
#Plot evolution of patients
plot_evolution <- function(df_plot, title){
  p <- ggplot(df_plot, aes(x = cycle, y = value, group = state, color = state)) +
    geom_line() + 
    ggtitle(title) +
    xlab("Months") + ylab("Probability") + 
    theme_minimal() +
    theme(
      plot.title = element_text(color="dodgerblue4", size=14, face="bold.italic", hjust = 0.5), # Center title
      axis.title.x = element_text(color="darkolivegreen4", size=14, face="bold"),
      axis.title.y = element_text(color="#993333", size=14, face="bold"),
      axis.line = element_line(color = "black"),  # Set color of axis lines
      legend.text = element_text(family = "Arial"),  # Set legend font family
      legend.title = element_text(family = "Arial", face = "bold")  # Set legend title font family
    ) +
    scale_x_continuous(breaks = seq(0, max(df_plot$cycle), by = 6)) +  # Adjust x-axis breaks
    scale_color_manual(values = c("Dead" = "red", "Progressed_drug" = "gold", "Progressed_nodrug" = "green", "ProgressedAE"="cyan", "ProgressedILD"= "blue", "ProgressionFree"= "purple", "ProgressionFreeAE"= "magenta", "ProgressionFreeILD"="pink"),
                       labels = c("Dead", "Progressed after drug", "Progressed after no drug", "Progressed after AE", "Progressed after ILD", "Progression-Free", "Progression-Free after AE", "Progression-Free after ILD"),
                       name = "State")  # Update legend title
  
  return(p)
}

# Clear plotting region
dev.off()

# Example calls to the function
plot5 <- plot_evolution(df_plot_tdxd_chemo, "Evolution of patients with T-Dxd first and chemo second")
plot6 <- plot_evolution(df_plot_chemo_chemo, "Evolution of patients with chemo first and second")
plot7 <- plot_evolution(df_plot_chemo_tdxd, "Evolution of patients with chemo first and T-Dxd second")
plot8 <- plot_evolution(df_plot_tdxd_sg, "Evolution of patients with T-Dxd first and SG second")

# Combine plots vertically
combined_plots <- arrangeGrob(plot6, plot5, ncol = 1)
print(combined_plots)
grid.draw(combined_plots)














#Plot efficiency frontiers
tdxd_icers <- calculate_icers(cost = df_res$DiscountedCost, 
                              effect = df_res$DiscountedQALY, 
                              strategies = df_res$Strategy)

plot(tdxd_icers) + 
  scale_x_continuous(n.breaks = 10) + 
  scale_y_continuous(n.breaks = 10) + 
  theme(legend.position = c(0.3, 0.7))




#One Way Sensitivity Analysis
one_way_sensitivity_tdxd_price <- function(df, dr_v){
  
  cost <- seq(1000, 15000, 500)
  
  n_sims = length(cost)
  
  icer <- c()
  
  tdxd_chemo_qaly <- c(qaly_pf = 0.65, 
                       qaly_p_drug = 0.54, 
                       qaly_p_nodrug = 0.26, 
                       qaly_pfAE = 0.547, 
                       qaly_pAE = 0.54, 
                       qaly_pfILD = 0.547, 
                       qaly_pILD = 0.54, 
                       decrement_qaly_ae = 0.05604)
  chemo_tdxd_qaly <- c(qaly_pf = 0.65, 
                       qaly_p_drug = 0.54, 
                       qaly_p_nodrug = 0.26, 
                       qaly_pfAE = 0.547, 
                       qaly_pAE = 0.54, 
                       qaly_pfILD = 0.547, 
                       qaly_pILD = 0.54, 
                       decrement_qaly_ae = 0.0714327)
  chemo_chemo_qaly <- c(qaly_pf = 0.65, 
                        qaly_p_drug = 0.54, 
                        qaly_p_nodrug = 0.26, 
                        qaly_pfAE = 0.547, 
                        qaly_pAE = 0.54, 
                        qaly_pfILD = 0.547, 
                        qaly_pILD = 0.54, 
                        decrement_qaly_ae = 0.0714327)
  
  df_tdxd_chemo = df[[1]]
  df_chemo_tdxd = df[[2]]
  
  # Initialize empty data frames to store results
  df_results <- data.frame()
  
  for(i in cost){
    
    tdxd_chemo_cost <- c(cost_pf = i, 
                         cost_p_drug = 7203.56, 
                         cost_p_nodrug = 10882.60, 
                         cost_pfAE = 5093.37, 
                         cost_pAE = 10882.60, 
                         cost_pfILD = 5093.37, 
                         cost_pILD = 10882.60,
                         additional_cost = 7519.15)
    
    chemo_tdxd_cost <- c(cost_pf = 7203.56, 
                         cost_p_drug = i, 
                         cost_p_nodrug = 10882.60, 
                         cost_pfAE = 5093.37, 
                         cost_pAE = 10882.60, 
                         cost_pfILD = 5093.37, 
                         cost_pILD = 10882.60,
                         additional_cost = 10435.36588)
  
      chemo_chemo_cost <- c(cost_pf = 7203.56, 
                          cost_p_drug = 5093.37, 
                          cost_p_nodrug = 10882.60, 
                          cost_pfAE = 5093.37, 
                          cost_pAE = 10882.60, 
                          cost_pfILD = 5093.37, 
                          cost_pILD = 10882.60,
                          additional_cost = 10435.36588)
    
    
      # Calculate summary data for each comparison
      res_vec_tdxd_chemo <- calc_summary_data(df = df_tdxd_chemo, 
                                              cost_data = tdxd_chemo_cost, 
                                              qaly_data = tdxd_chemo_qaly, 
                                              dr_v = dr_v)
      
      res_vec_chemo_tdxd <- calc_summary_data(df = df_chemo_tdxd, 
                                              cost_data = chemo_tdxd_cost, 
                                              qaly_data = chemo_tdxd_qaly, 
                                              dr_v = dr_v)
      
      res_vec_chemo_chemo <- calc_summary_data(df = df_chemo_tdxd, 
                                               cost_data = chemo_chemo_cost, 
                                               qaly_data = chemo_chemo_qaly, 
                                               dr_v = dr_v)
      
      # Calculate ICER for each scenario and store in the results vector
      icer_tdxd_chemo <- (res_vec_tdxd_chemo[2] - res_vec_chemo_tdxd[2]) / (res_vec_tdxd_chemo[4] - res_vec_chemo_tdxd[4])
      icer_tdxd_chemo_chemo <- (res_vec_tdxd_chemo[2] - res_vec_chemo_chemo[2]) / (res_vec_tdxd_chemo[4] - res_vec_chemo_chemo[4])
      
      icer <- c(icer, icer_tdxd_chemo, icer_tdxd_chemo_chemo)
      
      # Prepare data frame to store results
      df_results <- rbind(df_results,
                          data.frame(Willingness2Pay = c(icer_tdxd_chemo, icer_tdxd_chemo_chemo),
                                     Cost = rep(i, 2),
                                     Comparison = c('TDXd-> Chemo vs. Chemo->TDXd', 'TDXd-> Chemo vs. Chemo Only')))
      
  }
  
  # Plotting the results
  library(ggplot2)
  
  ggplot(data = df_results, aes(x = Willingness2Pay, y = Cost, color = Comparison, group = Comparison)) +
    geom_line() +
    geom_point() +
    labs(x = "ICER", y = "T-DXd Cost ($)", color = "Treatment Strategies Comparisons") +
    scale_x_continuous(limits = c(0, 300000), breaks = seq(0, 300000, 50000)) +
    scale_y_continuous(limits = c(5000, 12500), breaks = seq(5000, 12500, 2500)) +
    theme_minimal() +
    theme(axis.line = element_line(size = 1, color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = c(0.15, 0.85),
          legend.background = element_rect(fill = "transparent"),
          legend.box.margin = margin(5, 5, 5, 5),
          legend.title = element_text(family = "Arial", face = "bold"),
          legend.text = element_text(family = "Arial", face = "bold"),
          axis.text = element_text(family = "Arial", face = "bold"),
          axis.title = element_text(family = "Arial", face = "bold"))
}

#NEEDS WORK
# Example usage:
df_one_way = list(df_tdxd_chemo, df_chemo_tdxd, df_tdxd_sg)
one_way_sensitivity_tdxd_price(df_one_way, dr_v = df_list_tdxd_chemo[[3]])












source("PSA_8state.R")


#PSA plot
plot_psa_scatter <- function(df_psa_res, group_name1, group_name2, group_name3, group_name4){
  X <- split(df_psa_res, df_psa_res$group)
  means <- data.frame(
    mean_cost = c(mean(X[[group_name1]]["DiscountedCost"][[1]]), 
                  mean(X[[group_name2]]["DiscountedCost"][[1]]), 
                  mean(X[[group_name3]]["DiscountedCost"][[1]]), 
                  mean(X[[group_name4]]["DiscountedCost"][[1]])),
    mean_qaly = c(mean(X[[group_name1]]["DiscountedQALY"][[1]]), 
                  mean(X[[group_name2]]["DiscountedQALY"][[1]]), 
                  mean(X[[group_name3]]["DiscountedQALY"][[1]]), 
                  mean(X[[group_name4]]["DiscountedQALY"][[1]])),
    group = c(group_name1, group_name2, group_name3, group_name4)
  )
  
  # Define colors
  colors <- c("Chemo Only" = "red", 
              "Chemo->TDXd" = "deeppink",  
              "TDXd->Chemo" = "green", 
              "TDXd->SG" = "darkkhaki")
  
  # Define x and y axis limits
  x_limits <- c(130000, 400000)
  y_limits <- c(0, 1.3)
  
  plot2 <- ggplot(df_psa_res[df_psa_res$group == "Chemo Only" | df_psa_res$group == "TDXd->Chemo", ], 
                  aes(x = DiscountedCost, y = DiscountedQALY, col = group)) + 
    geom_point() + 
    geom_point(data = means[means$group == "Chemo Only" | means$group == "TDXd->Chemo", ],  
               mapping = aes(x = mean_cost, y = mean_qaly), size = 3, col = "black") +
    scale_color_manual(values = colors) +  
    labs(x = bquote(bold("Discounted Cost ($)")), y = bquote(bold("Discounted QALY"))) +
    ggtitle(bquote(bold("Probabilistic Sensitivity Analysis for Chemo Only vs. TDXd->Chemo"))) +
    theme_bw(base_size = 14) +
    theme(
      text = element_text(family = "Arial", face = "bold"),
      panel.grid = element_blank(),
      axis.line = element_line(),
      axis.text = element_text(size = 12, family = "Arial", face = "bold"),
      axis.title = element_text(size = 14, family = "Arial", face = "bold"),
      plot.title = element_text(size = 12, family = "Arial", face = "bold", hjust = 0.5),
      legend.position = "none"
    ) +
    scale_x_continuous(limits = x_limits) +
    scale_y_continuous(limits = y_limits, breaks = c(0.3, 0.6, 0.9, 1.2))
  
  plot3 <- ggplot(df_psa_res[df_psa_res$group == "Chemo->TDXd" | df_psa_res$group == "TDXd->Chemo", ], 
                  aes(x = DiscountedCost, y = DiscountedQALY, col = group)) + 
    geom_point() + 
    geom_point(data = means[means$group == "Chemo->TDXd" | means$group == "TDXd->Chemo", ],  
               mapping = aes(x = mean_cost, y = mean_qaly), size = 3, col = "black") +
    scale_color_manual(values = colors) +  
    labs(x = bquote(bold("Discounted Cost ($)")), y = bquote(bold("Discounted QALY"))) +
    ggtitle(bquote(bold("Probabilistic Sensitivity Analysis for Chemo->TDXd vs. TDXd->Chemo"))) +
    theme_bw(base_size = 14) +
    theme(
      text = element_text(family = "Arial", face = "bold"),
      panel.grid = element_blank(),
      axis.line = element_line(),
      axis.text = element_text(size = 12, family = "Arial", face = "bold"),
      axis.title = element_text(size = 14, family = "Arial", face = "bold"),
      plot.title = element_text(size = 12, family = "Arial", face = "bold", hjust = 0.5),
      legend.position = "none"
    ) +
    scale_x_continuous(limits = x_limits) +
    scale_y_continuous(limits = y_limits, breaks = c(0.3, 0.6, 0.9, 1.2))
  
  plot4 <- ggplot(df_psa_res[df_psa_res$group == "TDXd->Chemo" | df_psa_res$group == "TDXd->SG", ], 
                  aes(x = DiscountedCost, y = DiscountedQALY, col = group)) + 
    geom_point() + 
    geom_point(data = means[means$group == "TDXd->Chemo" | means$group == "TDXd->SG", ],  
               mapping = aes(x = mean_cost, y = mean_qaly), size = 3, col = "black") +
    scale_color_manual(values = colors) +  
    labs(x = bquote(bold("Discounted Cost ($)")), y = bquote(bold("Discounted QALY"))) +
    ggtitle(bquote(bold("Probabilistic Sensitivity Analysis for TDXd->Chemo vs. TDXd->SG"))) +
    theme_bw(base_size = 14) +
    theme(
      text = element_text(family = "Arial", face = "bold"),
      panel.grid = element_blank(),
      axis.line = element_line(),
      axis.text = element_text(size = 12, family = "Arial", face = "bold"),
      axis.title = element_text(size = 14, family = "Arial", face = "bold"),
      plot.title = element_text(size = 12, family = "Arial", face = "bold", hjust = 0.5),
      legend.position = "none"
    ) +
    scale_x_continuous(limits = x_limits) +
    scale_y_continuous(limits = y_limits, breaks = c(0.3, 0.6, 0.9, 1.2))
  
  # Remove unknown black dots if they represent the mean of all points
  plot2$layers <- plot2$layers[-c(2)]
  plot3$layers <- plot3$layers[-c(2)]
  plot4$layers <- plot4$layers[-c(2)]
  
  grid.arrange(plot2, plot3, plot4, nrow = 2, ncol = 2, heights = c(1, 1), widths = c(1, 1), layout_matrix = rbind(c(1, 2), c(3, NA)))
}

# Run the function
plot_psa_scatter(res[[1]], "Chemo Only", "Chemo->TDXd", "TDXd->Chemo", "TDXd->SG")











#CEAC plot
plot_ceac_curve <- function(df_psa_res, group_name1, group_name2, group_name3, group_name4){
  X<-split(df_psa_res, df_psa_res$group)
  Y = X$`TDxD-Chemo`
  
  
  x = seq(50000, 400000, 10000)
  y = c()
  for(i in x){
    y = append(y, sum(Y$'ICER'<i)/length(Y$'ICER'))
  }
  
  df <- data.frame(threshold = x,
                   probability = y)
  print(df)
  
  ggplot(data=df, aes(x=threshold, y=probability, group=1)) +
    geom_line()+
    geom_point()
}


plot_ceac_curve(res[[1]], "Chemo-Chemo","Chemo-TDxd", "TDxD-Chemo", "TDxD-SG")

