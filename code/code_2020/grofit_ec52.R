# this only pertains to ec52 run on 19 march 2020
# tobias mueller



#------------------- setup and packages ---------------------------------------

setwd("C:/Users/obiew/Desktop/vannette/96wellplate_growthassay")

library(grofit)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(data.table)



#---------------------- part I --  reformatting the df --------------------------------------------------
#                    

# read in data
#d <- read.csv("input/nectarchem_ec102_11jun2020_tm.csv")
#labels <- read.csv("input/labels.csv")

#read in data without control wells
d <- read.csv("input/ec52_19march2020_nocontrol.csv")
labels <- read.csv("input/labels_16mar2020_nocontrol.csv")

# fix labels classes
labels$treatment <- as.factor(labels$treatment)
labels$well <- as.factor(labels$well)
# get rid of replicate column - not sure why thats there.. *shake fist at past me* "arrrggg"
labels<- labels[,-3]
# make time in minutes rounded to nearest 15 
# (technically first read is at 13:28:00 but we'll call it 15 minutes for now)
d$time <- (seq.int(nrow(d)))*15

#now to get it to match the required df format

# transpose df
d <- t(d)
d <- as.data.frame(d)

# get the row labels as column 1
setDT(d, keep.rownames = TRUE)[]

# make time the headers
d <- as.data.frame(d)
# apparently it keeps de-dataframing into lists
names(d) <- d[1,]
d <- d[-1,]

# merge raw data from plate reader with treatment labels
gr <- merge(d, labels, by.x="time", by.y="well", incomparables=NA)

# then for my sanity we're gonna move the treatment labels to the front of the df
# not that it matters but I like to look at it from time to time

#now add 2 random columns because the package demands it
gr <- cbind(gr, microbe="ec52", date="march 19")

#useing function moveMe, move treatment information to the front

moveMe <- function(data, tomove, where = "last", ba = NULL) {
  temp <- setdiff(names(data), tomove)
  x <- switch(
    where,
    first = data[c(tomove, temp)],
    last = data[c(temp, tomove)],
    before = {
      if (is.null(ba)) stop("must specify ba column")
      if (length(ba) > 1) stop("ba must be a single character string")
      data[append(temp, values = tomove, after = (match(ba, temp)-1))]
    },
    after = {
      if (is.null(ba)) stop("must specify ba column")
      if (length(ba) > 1) stop("ba must be a single character string")
      data[append(temp, values = tomove, after = (match(ba, temp)))]
    })
  x
}

gr <- moveMe(gr, c("treatment", "microbe", "date"), "first")

write.csv(gr, "input/reformatted_ec52.csv", row.names = FALSE)

# okay we now have our data frame in the correct format to move forward
# to the next step
# which is more dataframe reformatting!!









# --------------------- begin  ------------------------------------------------------
# the below adapted from Rachel Vannettes grofit anaylsis


# check the data frame is the same - it should be from left to right:
# 3 treatment information columns at the start - there MUST be 3
# then a row labeled "time" that contains all the well numbers, each well number being a seperate row
# to the right of this the column headers are now the time stamps
# the cell information in each row is the OD reading for the well at that time stamp




gr<-read.csv("input/reformatted_ec52.csv", na.strings=c(".", "NA"))
t<-read.csv("input/reformatted_ec52.csv", header=F)
times<-t[1,5:ncol(t)]
times<-matrix(times, byrow=T, ncol=ncol(gr)-4, nrow=nrow(gr))
# this removes labels and makes a matrix of just time stamps
times<-data.frame(times)
# then dataframes it

# create a number string of ODinit; replicate into a matrix the size of your raw data
ODinit<-gr[,5]
ODmat<-matrix(ODinit, byrow=FALSE, ncol=ncol(gr)-4, nrow=nrow(gr))


#do an interesting transformation
#Divide values by ODinit add 1 and take the natural log
tOD<-log((gr[,5:ncol(gr)]/ODmat)+1)


#what happens if we dont transform
tOD2 <- gr[,5:ncol(gr)]-ODmat
#then any negatives set to zero
tOD2 <- replace(tOD2, tOD2 < 0, 0)

grow.m2<-cbind(gr[,1],gr[,2:3],tOD2)








# for now ive comemnted out the below - doesnt seem useful for the time being

# here we're subsetting for just the the control 
#grow.t<-grow.m[which(grow.m[,1]=="control with microbe"),]
# 
# then creating a matching length df
#times.t<-times[2:9,]


# set controls for how the gro.fit modeling runs
control1<-grofit.control(fit.opt="b", log.y.gc=FALSE, interactive=T)

# neg.nan.act       -- Logical, indicates wether the program should stop when negative growth values or NA values appear (TRUE). Otherwise the program removes this values silently (FALSE). Improper values may be caused by incorrect data or input errors. Default: FALSE.
# clean.bootstrap   -- Logical, determines if negative values which occur during bootstrap should be removed (TRUE) or kept (FALSE). Note: Infinite values are always removed. Default: TRUE.
# suppress.messages -- Logical, indicates wether grofit messages (information about current growth curve, EC50 values etc.) should be displayed (FALSE) or not (TRUE). This option is meant to speed up the processing of high throuput data. Note: warnings are still displayed. Default: FALSE.
# fit.opt           -- Indicates wether the program should perform a model fit ("m"), a spline fit ("s") or both ("b"). Default: "b".
# log.x.gc          -- Logical, indicates wether a ln(x+ 1) should be applied to the time data of the growth curves. Default: FALSE.
# log.y.gc          -- Logical, indicates wether a  ln(y+ 1)should be applied to the growth data of the growth curves. Default: FALSE.
# interactive       -- Logical, controls whether the fit of each growth curve is controlled manually by the user. Default: TRUE.
# nboot.gc          -- Number of bootstrap samples used for the model free growth curve fitting. Use nboot.gc=0 to disable the bootstrap. Default: 0.
# smooth.gc         -- Parameter describing the smoothness of the spline fit; usually (not necessary) in (0;1]. Set smooth.gc=NULL causes the program to query an optimal value via cross validation techniques. Note: This is partly experimental. In future improved implementations of the smooth.spline function may lead to different results. See documentation of the R function smooth.spline for further details. Especially for datasets with few data points the option NULL might result in a too small smoothing parameter, which produces an error in smooth.spline. In that case the usage of a fixed value is recommended. Default: NULL.
# model.type        -- Character vector giving the names of the parametric models which should be fitted to the data. Default: c("gompertz", "logistic", "gompertz.exp", "richards").
# have.atleast      -- Minimum number of different values for the response parameter one shoud have for estimating a dose response curve. Note: All fit procedures require at least six unique values. Default: 6.
# parameter         -- The column of the output table which should be used for creating a dose response curve. See documentation of gcFit, drFit or summary.gcFit for further details. Default: 9, which represents the maximum slope of the parametric growth curve fit.
# smooth.dr         -- Smoothing parameter used in the spline fit by smooth.spline during dose response curve estimation. Usually (not necessesary) in (0; 1]. See documentation of smooth.spline for further details. Default: NULL.
# log.x.dr          -- Logical, indicates wether a ln( x+1) should be applied to the concentration data of the dose response curves. Default: FALSE.
# log.y.dr          -- Logical, indicates wether a ln( y+ 1) should be applied to the response data of the dose response curves. Default: FALSE.
# nboot.dr          -- Numeric value, defining the number of bootstrap samples for EC50 estimation. Use nboot.dr=0 to disable bootstrapping. Default: 0.


#growth.test<-gcFit(times,grow.m, control=control1)
#the model fits are terrible with the log transform

growth.test<-gcFit(times,grow.m2, control=control1)
#if we dont transform and only subtract the starting values the fits are MUCH better 


parms<-summary.gcFit(growth.test)
#use the built in summary function and write to a df




# lambda = length of lag phase
# mu = maximum growth rate\
# A = carrying capacity



#then remove unreliable wells
parms<-parms[parms$reliability==TRUE,]

# #and set NAs to zero for nu and A
parms$mu.model[is.na(parms$mu.model)] <- 0
parms$A.model[is.na(parms$A.model)] <- 0

#####  adjusting to z value relative to control with microbe ------------------------------------

#ideally I should figure out these two movemes but they have slightly differnt syntax 
#so for now ill use both 
moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}







#create adjuste Z score column and move it next to model column for A, mu, lambda
parms$z.mu<-(parms$mu.model-mean(parms$mu.model[which(parms$TestId=="control with microbe")]))/sd(parms$mu.model[which(parms$TestId=="control with microbe")])
parms <- parms[moveme(names(parms), "z.mu before mu.model")]

parms$z.A<-(parms$A.model-mean(parms$A.model[which(parms$TestId=="control with microbe")]))/sd(parms$A.model[which(parms$TestId=="control with microbe")])
parms <- parms[moveme(names(parms), "z.A before A.model")]

parms$z.lambda<-(parms$lambda.model-mean(parms$lambda.model[which(parms$TestId=="control with microbe")]))/sd(parms$lambda.model[which(parms$TestId=="control with microbe")])
parms <- parms[moveme(names(parms), "z.lambda before lambda.model")]


# 
# parms$mu.adjusted <- parms$mu.model/mean(parms$mu.model[which(parms$TestId=="control with microbe")])
# parms <- parms[moveme(names(parms), "mu.adjusted before mu.model")]
# 
# #do this also for A
# parms$A.adjusted <- parms$A.model/mean(parms$A.model[which(parms$TestId=="control with microbe")])
# parms <- parms[moveme(names(parms), "A.adjusted before A.model")]
# 
# #and lambda
# parms$lambda.adjusted <- parms$lambda.model/mean(parms$lambda.model[which(parms$TestId=="control with microbe")])
# parms <- parms[moveme(names(parms), "lambda.adjusted before lambda.model")]

write.csv(parms, "outputs/ec52grofit_summary.csv")

# ---------------------- plots -------------------


#plot mu
p1 <- ggplot(parms, aes(x=TestId, y=mu.model)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(y=mu.model), color="red", size=.5) +
  theme(axis.text.x=element_text(angle = -70, hjust = 0))

p2 <- ggplot(parms, aes(x=TestId, y=z.mu)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(y=z.mu), color="red", size=.5) +
  theme(axis.text.x=element_text(angle = -70, hjust = 0))

#plot A
p3 <- ggplot(parms, aes(x=TestId, y=A.model)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(y=A.model), color="red", size=.5) +
  theme(axis.text.x=element_text(angle = -70, hjust = 0))

p4 <- ggplot(parms, aes(x=TestId, y=z.A)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(y=z.A), color="red", size=.5) +
  theme(axis.text.x=element_text(angle = -70, hjust = 0))

#plot lambda
p5 <- ggplot(parms, aes(x=TestId, y=lambda.model)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(y=lambda.model), color="red", size=.5) +
  theme(axis.text.x=element_text(angle = -70, hjust = 0))

p6 <- ggplot(parms, aes(x=TestId, y=z.lambda)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(y=z.lambda), color="red", size=.5) +
  theme(axis.text.x=element_text(angle = -70, hjust = 0))


windows()
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)

g <- arrangeGrob(p1, p2, p3, p4, p5, p6, nrow = 3)
ggsave("outputs/ec52plot.png", g, height = 8, width = 6, dpi = 600)







