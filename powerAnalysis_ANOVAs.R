setwd("C:/Users/Abrar/Dropbox/UNC_OneDrive/Shaikh Lab/Simulation Power Ananlysis/")
rm(list=ls())
library(pwr)
library(pwr2)
library(easypower)

#power analysis for 2-way ANOVA
ss.2way(a=4,b=4,alpha=0.05,beta=0.1,f.A = 0.4, f.B = 0.4, B=100)
#beta = 1- power - so when beta = 0.1 it means 90% power, and beta=0.2 means 80% power. 

#3-way ANOVA  
main.eff1 <- list(name = "Sex", levels = 2, eta.sq = "med")
main.eff2 <- list(name = "Time", levels = 3, eta.sq = "med")
main.eff3 <- list(name = "Treatment", levels = 4, eta.sq = "med")
n.multiway(iv1 = main.eff1, iv2 = main.eff2, iv3 = main.eff3, interaction.eta2 = 0.0588) #0.0588 = 0.06 which is assuming a moderate effect size for the interaction terms 
#Link that contains the 3 way ANOVA code: https://cran.r-project.org/web/packages/easypower/vignettes/User_Input.html 