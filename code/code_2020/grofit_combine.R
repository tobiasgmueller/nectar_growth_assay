#this script will take the outputs of each plates individual script and combine them
# tobias Mueller 8/2020


setwd("C:/Users/obiew/Desktop/vannette/96wellplate_growthassay")

# first read in packages
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(broom)
library(rstatix)
library(svglite)
#then read in data

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



#now lets lets drop down some graphs
# here the line at 1 is the average value for the control growth wells 

#for alpha
p1<-ggplot(df, aes(x=strain, y=z.A, color = strain )) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(y=z.A), color="red", size=.5) +
  facet_wrap(df$treatment, scales = "free")+
  geom_hline(aes(yintercept=0)) +
  theme(axis.text.x=element_text(angle = -70, hjust = 0))
 
ggsave("outputs/combined_A_zed.png", p1, height = 8, width = 10, dpi = 400)
ggsave("outputs/combined_A_zed.svg", p1, height = 8, width = 10, dpi = 400)


ggplot(df, aes(x=treatment, y=z.A, color = treatment )) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(y=z.A), color="red", size=.5) +
  facet_wrap(df$strain, scales = "free")+
  geom_hline(aes(yintercept=0)) +
  theme(axis.text.x=element_text(angle = -70, hjust = 0))


# ggplot(df, aes(x=treatment, y=A.adjusted, color = strain )) +
#   geom_boxplot(aes(), alpha=0.5) +
#   geom_point(aes(y=A.adjusted), color="red", size=.5) +
#   facet_wrap(df$strain, scales = "free")+
#   geom_hline(aes(yintercept=1)) +
#   theme(axis.text.x=element_text(angle = -70, hjust = 0))


# for mu
p2<- ggplot(df, aes(x=strain, y=z.mu, color = strain )) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(y=z.mu), color="red", size=.5) +
  facet_wrap(df$treatment, scales = "free")+
  geom_hline(aes(yintercept=0)) +
  theme(axis.text.x=element_text(angle = -70, hjust = 0))

ggsave("outputs/combined_mu_zed.png", p2, height = 8, width = 10, dpi = 400)
ggsave("outputs/combined_mu_zed.svg", p2, height = 8, width = 10, dpi = 400)
       


Pstrain <- ggplot(df, aes(x=treatment, y=mu.model, color = treatment )) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(y=mu.model), color="red", size=.5) +
  facet_wrap(df$strain, scales = "free")+
  geom_hline(aes(yintercept=0)) +
  theme(axis.text.x=element_text(angle = -70, hjust = 0))

ggsave("outputs/combined_mu_strain.svg", Pstrain, height = 8, width = 10, dpi = 400)

#for lambda
p3<-ggplot(df, aes(x=strain, y=z.lambda, color = strain )) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(y=z.lambda), color="red", size=.5) +
  facet_wrap(df$treatment, scales = "free")+
  geom_hline(aes(yintercept=0)) +
  theme(axis.text.x=element_text(angle = -70, hjust = 0))

ggsave("outputs/combined_lambda_zed.png", p3, height = 8, width = 10, dpi = 400)
ggsave("outputs/combined_lambda_zed.svg", p3, height = 8, width = 10, dpi = 400)


lambdastrain <- ggplot(df, aes(x=treatment, y=z.lambda, color = treatment )) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(y=z.lambda), color="red", size=.5) +
  facet_wrap(df$strain, scales = "free")+
  geom_hline(aes(yintercept=0)) +
  theme(axis.text.x=element_text(angle = -70, hjust = 0))

ggsave("outputs/combined_lambda_strain.svg", lambdastrain, height = 8, width = 10, dpi = 400)

windows()
grid.arrange(p1, p2, p3, nrow=1)



#anovas for mu
temp <- filter(df, strain == c("bb362"))
aovmubb362 <-aov(data=temp, mu.model ~ treatment)
TukeyHSD(aovmubb362)

temp <- filter(df, strain == c("ec102"))
aovmuec102 <-aov(data=temp, mu.model ~ treatment)
TukeyHSD(aovmuec102)

temp <- filter(df, strain == c("ec52"))
aovmuec52 <-aov(data=temp, mu.model ~ treatment)
TukeyHSD(aovmuec52)

temp <- filter(df, strain == c("ec84"))
aovmuec84 <-aov(data=temp, mu.model ~ treatment)
TukeyHSD(aovmuec84)

temp <- filter(df, strain == c("q1f2"))
aovmuq1f2 <-aov(data=temp, mu.model ~ treatment)
TukeyHSD(aovmuq1f2)


aovmu <-aov(data=temp, mu.model ~ treatment+strain)
TukeyHSD(aovmu)



# perhaps some simple t-tests to see what different from 0?
# the below subsets based on a combination of treament and strain and then runs t-tests across those groups
