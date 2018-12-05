setwd("C:/Users/Abrar/Dropbox/UNC_OneDrive/Shaikh Lab/Simulation Power Ananlysis/")
rm(list=ls())
library(pwr)

#determine effect sizes from previous study
data_t <- read.csv("Power_8.17.18.csv", header = TRUE) 

#calculating effect size for 5S-HETE/14-HDHA study using Cohen's d method
#d = mean(pre) - mean(post) / pooled SD 
#instead of pre and post it can just be mean 1 and mean 2
d_data_t <- (mean(data_t$Pre) - mean(data_t$Post[1:7])) / sqrt(((sd(data_t$Pre)^2 + sd(data_t$Post[1:7])^2)/2))

pwr.t.test(n = 70, d = 1.16, sig.level = 0.05, type = c("two.sample"))
#You can specify alternative="two.sided", "less", or "greater" to indicate a two-tailed, or one-tailed test. A two tailed test is the default

#Link with multiple power analyses: https://www.statmethods.net/stats/power.html 