

# this will be a quick test to see if decreases in growth curves correlate to stressful environments

setwd("C:/Users/obiew/Desktop/vannette/96wellplate_growthassay")

library(grofit)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(data.table)


# I think an easy way to do this would be to create a huge datafile with all the raw OD readings and the parametrics
# then for each well take the end measurement and subtract the peak measurement
# any negative values indicate a drop

# oof wait that probably wont work for ec84 where values just everywhere.... m
# maybe we take end measurements - A?


d1 <- read.csv("outputs/ec52grofit_summary.csv")
d2 <- read.csv("outputs/ec102grofit_summary.csv")
d3 <- read.csv("outputs/bb362grofit_summary.csv")
d4 <- read.csv("outputs/q1f2grofit_summary.csv")
d5 <- read.csv("outputs/ec84grofit_summary.csv")

#df should include all summary parametrics 
df <- Reduce(function(x, y) merge(x, y, all=TRUE), list(d1, d2, d3, d4, d5))
#then lets simplify these grofit column names

df <- rename(df, c("strain" = "AddId", "date" = "concentration", "treatment" = "TestId"))

df$treatment <- as.factor(df$treatment)
df$strain <- as.factor(df$strain)




raw1<-read.csv("input/reformatted_ec52.csv", na.strings=c(".", "NA"))
raw2<-read.csv("input/reformatted_ec102.csv", na.strings=c(".", "NA"))
raw3<-read.csv("input/reformatted_bb362.csv", na.strings=c(".", "NA"))
raw4<-read.csv("input/reformatted_ec84.csv", na.strings=c(".", "NA"))
raw5<-read.csv("input/reformatted_q1f2.csv", na.strings=c(".", "NA"))

df2 <- Reduce(function(x, y) merge(x, y, all=TRUE), list(raw1, raw2, raw3, raw4, raw5))



dfall <- join(df, df2, by= c("treatment", "strain" = "microbe",    ), type = "left", match = "first")