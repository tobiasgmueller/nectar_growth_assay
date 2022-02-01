## this is a test of the growthrates package for modeling growth curves

## 6/25/20 Tobias Mueller


# ill set up the same wd as my other work and pull those input files for this
setwd("C:/Users/obiew/Desktop/vannette/96wellplate_growthassay")

library(growthrates)
library(lattice)


library(growthcurver)
library(lubridate)
library(magrittr)
library(reshape2)
library(ggplot2)
library(readr)
library(dplyr)
library(broom)
library(emmeans)


#lets use the rerun of ec102 A. pallulans - june 11th
d <- read.csv("input/nectarchem_ec102_11jun2020_tm.csv")
d_nocont <- read.csv("input/nectarchem_ec102_11jun2020_tm_nocontrol.csv")
labels <- read.csv("input/labels.csv")
labels_no_control <- read.csv("input/labels_nocontrol.csv")






#so it appears to use growthrates you need a long dataframe so first we meltify

meltd <- melt(d, id.vars = "time", variable.name = "sample")

#now lets add our labels (also required I beleive)
df <- merge(meltd, labels, by.x="sample", by.y="well", incomparables=NA)
df$treatment <- as.factor(df$treatment)

#okay so it seems we need a running increasing integer for time with no duplicates
df$tim <- sequence(rle(as.character(df$sample))$lengths)

#lets try to use their recommendation of lattice for graphing
xyplot(value ~ tim|treatment, data = df,
       groups = replicate, pch = 16, cex = 0.5)

#Whoa... those are some wonky looking graphs.. but looks like it worked

# damn do these graphs take FOREVER to render









#### first ill try the "fitting models to single data sets" using the "easy linear method"
#which I think is one single curve


#first we spit out data into a list where each iten in the list is a single well
split.data <- multisplit(df, c("treatment", "sample"))
#then write a well to a df (in this case item #3 from list which equals well A3) 
dat <- split.data[[6]]

#then we fit
fit <- fit_easylinear(dat$tim, dat$value)
#im getting an "essentially perfect fit" warning... which is uh bad?

plot(fit)
#lol well im confused where this "essentially perfect" is coming from...

#okay so that plots a straight linear fit







#### now to step 2 - learning their non parametric fitting


#my understsanding from reading is I pick a starting value for p, mu, and k
#and then optionally i can pick limits for them (example of limits below)
p     <- c(y0 = .07, mumax = .2, K = .23)

#lower <- c(y0 = 1e-5, mumax = 0,   K = 0)
#upper <- c(y0 = .1, mumax = .6,   K = .6)

fit1 <- fit_growthmodel(FUN = grow_logistic, p = p, dat$tim, dat$value)

#fit1 <- fit_growthmodel(FUN = grow_logistic, p = p, dat$tim, dat$value, lower = lower, upper = upper)

plot(fit1)

# hmm.. it does plot a "fit" of sorts.. which changes drastically based on starting p values








#ahhh wait wait wait --  I found the "growth_twostep" - which is I think what I want 


p     <- c(yi = 0.07, ya = 0.052, kw = 0.1, mumax = 0.2, K = 0.23)
lower <- c(yi = 0, ya = 0, kw = 0,    mumax = 0,   K = 0)
upper <- c(yi = .1, ya = .5, kw = 1,   mumax = 3,   K = .3)

fit2 <- fit_growthmodel(FUN = grow_twostep, p = p, time = dat$tim, y = dat$value,
                        lower = lower, upper = upper)

plot(fit2)
# much nicer. But what happens when I change my parameters.. 
#ah everything changes
#I think still I want a package to determine those parameters for me... otherwise work 


#consensus on package: meh. Maybe if I was better at math?
