#This code was created by Marcus Moen, MS and edited by James Dickerson, MD, MS. and ***
#Oversight from Jeremy Goldhaber-Feibert, PhD and Fernando Alarid-Escudero, Ph.D.

options(scipen=999)
setwd("/Users/jamesdickerson/Library/CloudStorage/Box-Box/Dickerson Lab/Dickerson_Lab_Github/CEA_of_TDXd/")

#Plot file
library(ggplot2)
library(ggpubr)
library(dampack)













source("finding_transition_probabilities_7state.R")

source("finding_transition_probabilities_8state.R")

#Store the transition matrices from the "finding_transition_probabilities_8state.R" file
A_tdxd_chemo <- tm_tdxd_chemo
A_chemo_chemo <- tm_chemo_chemo
A_chemo_tdxd <- tm_chemo_tdxd
A_tdxd_sg <- tm_tdxd_sg


###TODO!!! Add the two lines for the two others as well. Combine the two plots on the first row with these.
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
  
  
  n = max(length(km_os_chemo_chemo), length(km_pf_chemo_chemo), length(km_os_tdxd_chemo), length(km_pf_tdxd_chemo))
  
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
  
  ae_test_chemo <- c(0)
  ae_test_tdxd <- c(0)
  ae_test_ild <- c(0)
  
  
  for(t in 1:n){

    s1_chemo_chemo <- s0_chemo_chemo %*% A_chemo_chemo[[t]]
    s1_tdxd_chemo <- s0_tdxd_chemo %*% A_tdxd_chemo[[t]]
    s1_chemo_tdxd <- s0_chemo_tdxd %*% A_chemo_tdxd[[t]]
    s1_tdxd_sg <- s0_tdxd_sg %*% A_tdxd_sg[[t]]
    
    #survival_chemo <- s1_chemo[1]+s1_chemo[4]
    #km_os_chemo_model <- append(km_os_chemo_model, survival_chemo)
    #km_pf_chemo_model <- append(km_pf_chemo_model, s1_chemo[1])
    
    #survival_tdxd <- s1_tdxd[1]+s1_tdxd[4]
    #km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
    #km_pf_tdxd_model <- append(km_pf_tdxd_model, s1_tdxd[1])
    
    survival_chemo_chemo <- km_os_chemo_chemo_model[t]*(s1_chemo_chemo[1]+s1_chemo_chemo[4]+s1_chemo_chemo[5])/(s1_chemo_chemo[1]+s1_chemo_chemo[4]+s1_chemo_chemo[5]+s0_chemo_chemo[1]*A_chemo_chemo[[t]][1,8]+s0_chemo_chemo[4]*A_chemo_chemo[[t]][4,8]+s0_chemo_chemo[5]*A_chemo_chemo[[t]][5,8])
    km_os_chemo_chemo_model <- append(km_os_chemo_chemo_model, survival_chemo_chemo)
    km_pf_chemo_chemo_model <- append(km_pf_chemo_chemo_model, km_pf_chemo_chemo_model[t]*s1_chemo_chemo[1]/(s1_chemo_chemo[1]+s0_chemo_chemo[1]*A_chemo_chemo[[t]][1,4]+s0_chemo_chemo[1]*A_chemo_chemo[[t]][1,5]+s0_chemo_chemo[1]*A_chemo_chemo[[t]][1,8]))
    
    survival_tdxd_chemo <- km_os_tdxd_chemo_model[t]*(s1_tdxd_chemo[1]+s1_tdxd_chemo[4]+s1_tdxd_chemo[5])/(s1_tdxd_chemo[1]+s1_tdxd_chemo[4]+s1_tdxd_chemo[5]+s0_tdxd_chemo[1]*A_tdxd_chemo[[t]][1,8]+s0_tdxd_chemo[4]*A_tdxd_chemo[[t]][4,8]+s0_tdxd_chemo[5]*A_tdxd_chemo[[t]][5,8])
    km_os_tdxd_chemo_model <- append(km_os_tdxd_chemo_model, survival_tdxd_chemo)
    km_pf_tdxd_chemo_model <- append(km_pf_tdxd_chemo_model, km_pf_tdxd_chemo_model[t]*s1_tdxd_chemo[1]/(s1_tdxd_chemo[1]+s0_tdxd_chemo[1]*A_tdxd_chemo[[t]][1,4]+s0_tdxd_chemo[1]*A_tdxd_chemo[[t]][1,5]+s0_tdxd_chemo[1]*A_tdxd_chemo[[t]][1,8]))
    
    survival_chemo_tdxd <- km_os_chemo_tdxd_model[t]*(s1_chemo_tdxd[1]+s1_chemo_tdxd[4]+s1_chemo_tdxd[5])/(s1_chemo_tdxd[1]+s1_chemo_tdxd[4]+s1_chemo_tdxd[5]+s0_chemo_tdxd[1]*A_chemo_tdxd[[t]][1,8]+s0_chemo_tdxd[4]*A_chemo_tdxd[[t]][4,8]+s0_chemo_tdxd[5]*A_chemo_tdxd[[t]][5,8])
    km_os_chemo_tdxd_model <- append(km_os_chemo_tdxd_model, survival_chemo_tdxd)
    km_pf_chemo_tdxd_model <- append(km_pf_chemo_tdxd_model, km_pf_chemo_tdxd_model[t]*s1_chemo_tdxd[1]/(s1_chemo_tdxd[1]+s0_chemo_tdxd[1]*A_chemo_tdxd[[t]][1,4]+s0_chemo_tdxd[1]*A_chemo_tdxd[[t]][1,5]+s0_chemo_tdxd[1]*A_chemo_tdxd[[t]][1,8]))
    
    survival_tdxd_sg <- km_os_tdxd_sg_model[t]*(s1_tdxd_sg[1]+s1_tdxd_sg[4]+s1_tdxd_sg[5])/(s1_tdxd_sg[1]+s1_tdxd_sg[4]+s1_tdxd_sg[5]+s0_tdxd_sg[1]*A_tdxd_sg[[t]][1,8]+s0_tdxd_sg[4]*A_tdxd_sg[[t]][4,8]+s0_tdxd_sg[5]*A_tdxd_sg[[t]][5,8])
    km_os_tdxd_sg_model <- append(km_os_tdxd_sg_model, survival_tdxd_sg)
    km_pf_tdxd_sg_model <- append(km_pf_tdxd_sg_model, km_pf_tdxd_sg_model[t]*s1_tdxd_sg[1]/(s1_tdxd_sg[1]+s0_tdxd_sg[1]*A_tdxd_sg[[t]][1,4]+s0_tdxd_sg[1]*A_tdxd_sg[[t]][1,5]+s0_tdxd_sg[1]*A_tdxd_sg[[t]][1,8]))
    
    s0_chemo_chemo <- s1_chemo_chemo
    s0_tdxd_chemo <- s1_tdxd_chemo
    s0_chemo_tdxd <- s1_chemo_tdxd
    s0_tdxd_sg <- s1_tdxd_sg
    
  }
  m = min(length(km_os_chemo_chemo), length(km_pf_chemo_chemo), length(km_os_tdxd_chemo), length(km_pf_tdxd_chemo))
  idx <- 1:m
  df <- data.frame(idx, km_os_chemo_chemo_model[1:m], km_os_tdxd_chemo_model[1:m], km_os_chemo_chemo[1:m])
  
  
  plot1 <- ggplot(df, aes(x=idx)) + 
    geom_line(aes(y = km_os_chemo_chemo_model[1:m]), color = "red") + 
    geom_line(aes(y = km_os_chemo_chemo[1:m]), color="blue", linetype="twodash") +
    geom_line(aes(y = km_os_tdxd_chemo_model[1:m]), color = "chocolate") + 
    geom_line(aes(y = km_os_tdxd_chemo[1:m]), color="green", linetype="twodash") +
    geom_line(aes(y = km_os_chemo_tdxd_model[1:m]), color = "deeppink") + 
    geom_line(aes(y = km_os_tdxd_sg_model[1:m]), color = "darkkhaki") + 
    ylim(0, 1) +
    ylab("Probability") +
    xlab("Month") +
    ggtitle("OS")+
    scale_x_continuous(breaks = round(seq(0, 18, by = 3),1))+theme(axis.text=element_text(size=12),
                                                                   axis.title=element_text(size=14))
  
  
  plot2 <- ggplot(df, aes(x=idx)) + 
    geom_line(aes(y = km_pf_chemo_chemo_model[1:m]), color = "red") + 
    geom_line(aes(y = km_pf_chemo_chemo[1:m]), color="blue", linetype="twodash") +
    geom_line(aes(y = km_pf_tdxd_chemo_model[1:m]), color = "chocolate") + 
    geom_line(aes(y = km_pf_tdxd_chemo[1:m]), color="green", linetype="twodash") +
    geom_line(aes(y = km_pf_chemo_tdxd_model[1:m]), color = "deeppink") + 
    geom_line(aes(y = km_pf_tdxd_sg_model[1:m]), color = "darkkhaki") + 
    ylim(0, 1) +
    ylab("Probability") +
    xlab("Month") +
    ggtitle("PF")+
    scale_x_continuous(breaks = round(seq(0, 18, by = 3),1))+theme(axis.text=element_text(size=12),
                                                                   axis.title=element_text(size=14))
  
 
  print(plot1)
  
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
  p <- ggplot(df_plot, aes(x = cycle,y = value, group = state, color = state)) +
    geom_line() + 
    ggtitle(title) +
    xlab("Months from start") + ylab("Probability") + 
    theme(
      plot.title = element_text(color="dodgerblue4", size=14, face="bold.italic"),
      axis.title.x = element_text(color="darkolivegreen4", size=14, face="bold"),
      axis.title.y = element_text(color="#993333", size=14, face="bold")
    )
  return(p)
}

plot_evolution(df_plot_tdxd_chemo, "Evolution of patients with TDxD first and chemo second")
plot_evolution(df_plot_chemo_chemo, "Evolution of patients with chemo first and second")
plot_evolution(df_plot_chemo_tdxd, "Evolution of patients with chemo first and TDxD second")
plot_evolution(df_plot_tdxd_sg, "Evolution of patients with TDxD first and SG second")



#Plot efficiency frontiers
tdxd_icers <- calculate_icers(cost = df_res$DiscountedCost, 
                              effect = df_res$DiscountedQALY, 
                              strategies = df_res$Strategy)

plot(tdxd_icers) + 
  scale_x_continuous(n.breaks = 10) + 
  scale_y_continuous(n.breaks = 10) + 
  theme(legend.position = c(0.3, 0.7))






source("PSA_8state.R")


#PSA plot
plot_psa_scatter <- function(df_psa_res, group_name1, group_name2, group_name3, group_name4){
  X<-split(df_psa_res, df_psa_res$group)
  means <- data.frame(mean_cost = c(mean(X[[group_name1]]["DiscountedCost"][[1]]), mean(X[[group_name2]]["DiscountedCost"][[1]]), mean(X[[group_name3]]["DiscountedCost"][[1]]), mean(X[[group_name4]]["DiscountedCost"][[1]])),
                      mean_qaly = c(mean(X[[group_name1]]["DiscountedQALY"][[1]]), mean(X[[group_name2]]["DiscountedQALY"][[1]]), mean(X[[group_name3]]["DiscountedQALY"][[1]]), mean(X[[group_name4]]["DiscountedQALY"][[1]])),
                      group = c(group_name1, group_name2, group_name3, group_name4))
  
  
  plot1 <- ggplot(df_psa_res, aes(x=DiscountedCost, y=DiscountedQALY, col = group)) + 
    geom_point() + 
    geom_point(data = means,  
               mapping = aes(x = mean_cost, y = mean_qaly, size = 1)) +
    labs(x = "Discounted Cost", y = "Discounted QALY") +
    theme_bw(base_size = 14)
  
  plot2 <- ggplot(df_psa_res[df_psa_res$group == "Chemo-Chemo"|df_psa_res$group=="TDxD-Chemo",], aes(x=DiscountedCost, y=DiscountedQALY, col = group)) + 
    geom_point() + 
    geom_point(data = means[means$group == "Chemo-Chemo"|means$group=="TDxD-Chemo",],  
               mapping = aes(x = mean_cost, y = mean_qaly, size = 1)) +
    labs(x = "Discounted Cost", y = "Discounted QALY") +
    theme_bw(base_size = 14)
  
  plot3 <- ggplot(df_psa_res[df_psa_res$group == "Chemo-TDxd"|df_psa_res$group=="TDxD-Chemo",], aes(x=DiscountedCost, y=DiscountedQALY, col = group)) + 
    geom_point() + 
    geom_point(data = means[means$group == "Chemo-TDxd"|means$group=="TDxD-Chemo",],  
               mapping = aes(x = mean_cost, y = mean_qaly, size = 1)) +
    labs(x = "Discounted Cost", y = "Discounted QALY") +
    theme_bw(base_size = 14)
  
  plot4 <- ggplot(df_psa_res[df_psa_res$group == "TDxD-Chemo"|df_psa_res$group=="TDxD-SG",], aes(x=DiscountedCost, y=DiscountedQALY, col = group)) + 
    geom_point() + 
    geom_point(data = means[means$group == "TDxD-Chemo"|means$group=="TDxD-SG",],  
               mapping = aes(x = mean_cost, y = mean_qaly, size = 1)) +
    labs(x = "Discounted Cost", y = "Discounted QALY") +
    theme_bw(base_size = 14)
  
  print(plot1)
  print(plot2)
  print(plot3)
  print(plot4)
}


plot_psa_scatter(res[[1]], "Chemo-Chemo","Chemo-TDxd", "TDxD-Chemo", "TDxD-SG")



#CEAC plot
plot_ceac_curve <- function(df_psa_res, group_name1, group_name2, group_name3, group_name4){
  X<-split(df_psa_res, df_psa_res$group)
  Y = X$`TDxD-Chemo`
  
  
  x = seq(50000, 400000, 10000)
  y = c()
  for(i in x){
    y = append(y, sum(Y$'ICER'<i)/length(Y$'ICER'))
  }
  
  df <- data.frame(treashold = x,
                   probability = y)
  print(df)
  
  ggplot(data=df, aes(x=treashold, y=probability, group=1)) +
    geom_line()+
    geom_point()
  
  
}


plot_ceac_curve(res[[1]], "Chemo-Chemo","Chemo-TDxd", "TDxD-Chemo", "TDxD-SG")


