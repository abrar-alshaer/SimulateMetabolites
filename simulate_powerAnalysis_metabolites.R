#Abrar Al-Shaer, Raz Shaikh Lab
#8/12/2018
#This program was adapted from Dr. Anthony Fodor, UNCC Bioinformatics Dept.

# this power simulation assumes perfectly a normal distribution 
# and then applies a t-test.
# Real data can have a different distribution
# so these results should be considered a rough approximation

setwd("C:/Users/Abrar/OneDrive - University of North Carolina at Chapel Hill/Shaikh Lab/Simulation Power Ananlysis/")
rm(list=ls())

#determine effect sizes from previous study
hete <- read.csv("5S-HETE study.csv", header = TRUE) #5S-HETE Study
hdha <- read.csv("14-HDHA study.csv", header = TRUE) #14-HDHA Study

#calculating effect size for 5S-HETE/14-HDHA study using Cohen's d method
#d = mean(pre) - mean(post) / pooled SD
d_hete <- (mean(hete$X5S.HETE_Pre) - mean(hete$X5S.HETE_Post)) / sqrt(((sd(hete$X5S.HETE_Pre)^2 + sd(hete$X5S.HETE_Post)^2)/2))
d_hdha <- (mean(hdha$X14.HDHA_Pre) - mean(hdha$X14.HDHA_Post)) / sqrt(((sd(hdha$X14.HDHA_Pre)^2 + sd(hdha$X14.HDHA_Post)^2)/2))

sampleSizeForEachGroup <- 100
numHypotheses <- 80 #same as # of metabolites
numSimulations <- 1000
fractionTruePositives <- 0.25
effectSize <- 0.54 #modeling for 5s-hete
#limitation of this simulations is that it only models based on the effect size of 14-hdha and 5s-hete, whereas we will 
#be modeling many different metabolites with varying effect sizes. However we can compensate with that by only specifying a 
#significant difference of 25% of our metabolites in the simulation and 75% as null, therefore increasing our power and needed sample size.

pValues <- vector();
truePositives <- vector()
power_bonf_ten <- vector()
power_bonf_five <- vector()
power_bh_ten <- vector()
power_bh_five <- vector()

for (k in 1:numSimulations) {
  numUnder <- 0
  numBonferroniUnder <-0
  numBonf <- 0
  numTruePositves <-0
  bhUnder <- 0
  benjHochUnder <- 0
for( j in 1 : numHypotheses ) 
{
  isATruePositive <- ( runif(1) <= fractionTruePositives )
  
  if( isATruePositive ) 	
    numTruePositves= numTruePositves + 1
  
  data<-vector()
  data2<-vector()
  
  for( i in 1 : sampleSizeForEachGroup)
  {
    data[i] = rnorm(1) #add regardless if it's a true positive
    
    if( isATruePositive  ) #add to data2 ONLY if it is a true positive
    {
      data2[i] = rnorm(1, mean=effectSize)
    }
    else
    {
      data2[i] = rnorm(1)
    }
  }
  
  pValues[j] <- t.test(data,data2)$p.value #is there is a significant difference between data and data2
  truePositives[j] <- isATruePositive
  
  if( isATruePositive  & pValues[j] < 0.05 )
    numUnder = numUnder + 1
  
  # if( isATruePositive  & pValues[j] < 0.05 / numHypotheses )
  #   numBonferroniUnder = numBonferroniUnder + 1	
  
} 

numTruePositves
numUnder/numTruePositves

myFrame <- data.frame(pValues,truePositives  )
myFrame<- myFrame[order(myFrame$pValues),]
myFrame$pAdjust <- p.adjust(myFrame$pValues, method="BH")
myFrame$Bonf <- p.adjust(myFrame$pValues, method="bonferroni")

for( a in 1:nrow(myFrame))
{
  if( myFrame$truePositives[a] & myFrame$pAdjust[a] < 0.10  )
    bhUnder <- bhUnder + 1
  
  if( myFrame$truePositives[a] & myFrame$pAdjust[a] < 0.05  )
    benjHochUnder <- benjHochUnder + 1
  
  if( myFrame$truePositives[a] & myFrame$Bonf[a] < 0.10  )
    numBonferroniUnder = numBonferroniUnder + 1
  
  if( myFrame$truePositives[a] & myFrame$Bonf[a] < 0.05  )
    numBonf = numBonf + 1
} 

power_bonf_ten[k] <- numBonferroniUnder/numTruePositves
power_bonf_five[k] <- numBonf/numTruePositves
#print(numBonferroniUnder/numTruePositves)

power_bh_ten[k] <- bhUnder / numTruePositves
power_bh_five[k] <- benjHochUnder / numTruePositves
#print(bhUnder / numTruePositves) 
}

mean(power_bonf_ten)
mean(power_bonf_five)
mean(power_bh_ten)
mean(power_bh_five)

