# Analysis of microbial growth in nectar with different compounds and concentrations

#### setup ######
setwd("C:/Users/obiew/Desktop/vannette/96wellplate_growthassay")

# read in data
gr<-read.csv("input/rlv_microbe_gr1_2.csv", na.strings=c(".", "NA"))
t<-read.csv("input/rlv_microbe_gr1_2.csv", header=F)
times<-t[1,5:262]
times<-matrix(times, byrow=T, ncol=ncol(gr)-4, nrow=nrow(gr))
times<-data.frame(times)
#install.packages("grofit")
library(grofit)
library(ggplot2)
# library(gdata)
# library(sciplot)
##################


# following previous analyses

##create a column for ODinit; replicate into a matrix
ODinit<-gr[,5]
ODmat<-matrix(ODinit, byrow=FALSE, ncol=ncol(gr)-4, nrow=nrow(gr))


#Divide values by ODinit and take the natural log
tOD<-log((gr[,5:ncol(gr)]/ODmat)+1)
grow.m<-cbind(gr[,1],gr[,2:3],tOD)
#here we're subsetting for just this microbe
grow.t<-grow.m[which(grow.m[,1]=="P. fluorescens"),]
#then creating a matching length df
times.t<-times[5:34,]


# run the growth models
yeast.control<-grofit.control(fit.opt="b", log.y.gc=FALSE, interactive=T)



par(mfrow=c(2,3))
growth.test<-gcFit(times,grow.m, control=yeast.control)


parms<-summary.gcFit(growth.test)



write.csv(parms, "TN_growth1_summary.csv")

growth.test<-gcFit(times.t,grow.t)



parms<-read.csv("TN_growth1_summary.csv")

#make new df with only reliable wells
parms.g1<-parms[parms$reliability==TRUE,]

#in parms write a 0 for mu on non-reliable wells
parms$mu.model[which(parms$reliability==FALSE)]<-0

#create the same df with only reliable wells
parms.g<-parms[parms$reliability==TRUE,]


#parms.g<-parms.g[parms.g$TestId!="Control",]
#parms$TestId<-drop.levels(parms$TestId)



#Examine growth rate
library(sciplot)
parms.nop<-parms[parms$TestId!="P. fluorescens",]
library(gdata)
parms.nop$TestId<-drop.levels(parms.nop$TestId)

windows()
par(mfrow=c(3,2))
#par has got to be the coolest function ever!! - my mind is blown


lineplot.CI(concentration, mu.model,group=TestId, data=parms.nop, subset=AddId=="Aucubin"|AddId=="Control", 
            legend=T, main="Aucubin", ylab="Max. growth rate (mu)", 
            xlab="Concentration")
lineplot.CI(concentration, mu.model,group=TestId, data=parms.nop, subset=AddId=="Catalpol"|AddId=="Control", 
            legend=T, main="Catalpol", ylab="Max. growth rate (mu)", 
            xlab="Concentration")
lineplot.CI(concentration, mu.model,group=TestId, data=parms.nop, subset=AddId=="Caffeine"|AddId=="Control", 
            legend=T, main="Caffeine", ylab="Max. growth rate (mu)", 
            xlab="Concentration")
lineplot.CI(concentration, mu.model,group=TestId, data=parms.nop, subset=AddId=="Nicotine"|AddId=="Control", 
            legend=T, main="Nicotine", ylab="Max. growth rate (mu)", 
            xlab="Concentration")
lineplot.CI(concentration, mu.model,group=TestId, data=parms.nop, subset=AddId=="Ouabain"|AddId=="Control", 
            legend=T, main="Ouabain", ylab="Max. growth rate (mu)", 
            xlab="Concentration")


par(mfrow=c(5,1), mar=c(3,2,1,0))
bargraph.CI(TestId, mu.model,group=concentration, data=parms.nop, subset=AddId=="Caffeine"|AddId=="Control", 
            legend=T , main="Caffeine")
bargraph.CI(TestId, mu.model,group=concentration, data=parms.nop, subset=AddId=="Nicotine"|AddId=="Control", 
            legend=T, main="Nicotine")
bargraph.CI(TestId, mu.model,group=concentration, data=parms.nop, subset=AddId=="Aucubin"|AddId=="Control", 
            legend=T , main="Aucubin")
bargraph.CI(TestId, mu.model,group=concentration, data=parms.nop, subset=AddId=="Catalpol"|AddId=="Control", 
            legend=T , main="Catalpol")
bargraph.CI(TestId, mu.model,group=concentration, data=parms.nop, subset=AddId=="Ouabain"|AddId=="Control", 
            legend=T , main="Ouabain")





#Examine mu (best-fit model) 
quartz()
par(mfrow=c(3,1))
bargraph.CI(AddId, mu.model,group=concentration, data=parms, subset=TestId=="C. rancensis", 
            legend=T, main="C. rancensis")
bargraph.CI(AddId, mu.model,group=concentration, data=parms, subset=TestId=="M. reukaufii", 
            legend=T, main="M. reukaufii")

bargraph.CI(AddId, mu.model,group=concentration, data=parms, subset=TestId=="P. fluorescens", 
            legend=T, main="P. fluorescens" )

#Examine A (best-fit model) 

quartz()
par(mfrow=c(3,2))
lineplot.CI(concentration, log(A.model),group=TestId, data=parms, subset=AddId=="Aucubin"|AddId=="Control", 
            legend=T, main="Aucubin", ylab="Max. OD (A)", 
            xlab="Concentration")
lineplot.CI(concentration, log(A.model),group=TestId, data=parms, subset=AddId=="Catalpol"|AddId=="Control", 
            legend=T, main="Catalpol", ylab="Max. OD (A)", 
            xlab="Concentration")
lineplot.CI(concentration, log(A.model),group=TestId, data=parms, subset=AddId=="Caffeine"|AddId=="Control", 
            legend=T, main="Caffeine", ylab="Max. OD (A)", 
            xlab="Concentration")
lineplot.CI(concentration, log(A.model),group=TestId, data=parms, subset=AddId=="Nicotine"|AddId=="Control", 
            legend=T, main="Nicotine", ylab="Max. OD (A)", 
            xlab="Concentration")
lineplot.CI(concentration, log(A.model),group=TestId, data=parms, subset=AddId=="Ouabain"|AddId=="Control", 
            legend=T, main="Ouabain", ylab="Max. OD (A)", 
            xlab="Concentration")

# Statistical analyses for maximum OD
fit<-lm(A.model~TestId*concentration, data=parms.g, subset=parms.g$AddId=="Ouabain"|parms.g$AddId=="Control")
summary(fit)
drop1(fit, .~., test="F")
fit<-lm(A.model~concentration, data=parms.g, subset=parms.g$AddId=="Catalpol"|parms.g$AddId=="Control")
summary(fit)
drop1(fit, .~., test="F")
fit<-lm(A.model~concentration, data=parms.g, subset=parms.g$AddId=="Nicotine"|parms.g$AddId=="Control")
summary(fit)
drop1(fit, .~., test="F")
fit<-lm(A.model~TestId+concentration, data=parms.g, subset=parms.g$AddId=="Caffeine"|parms.g$AddId=="Control")
summary(fit)
drop1(fit, .~., test="F")
fit<-lm(A.model~TestId*concentration, data=parms.g, subset=parms.g$AddId=="Aucubin"|parms.g$AddId=="Control")
summary(fit)
drop1(fit, .~., test="F")

# Statistical analyses for maximum growth rate
# all together

fit<-lm(mu.model~TestId*AddId*(concentration), data=parms.nop)
summary(fit)
drop1(fit, .~., test="F")
fit<-lm(A.model~TestId*AddId+(concentration), data=parms.nop)
fit<-lm(mu.model~TestId*AddId*(concentration), data=parms.nop, 
        subset=concentration!=5)



quartz()
lineplot.CI(concentration, mu.model, TestId, data=parms.nop)
lineplot.CI(concentration, A.model, TestId, data=parms.nop)


# separately


fit<-lm(mu.model~TestId+(concentration), data=parms.nop, 
        subset=parms.nop$AddId=="Ouabain"|parms.nop$AddId=="Control")
summary(fit)
drop1(fit, .~., test="F")

# No effect of ouabain


fit<-lm(mu.model~TestId+(concentration), data=parms.nop, 
        subset=parms.nop$AddId=="Aucubin"|parms.nop$AddId=="Control")
summary(fit)
drop1(fit, .~., test="F")
#interaction is significant: effect on P. fluorescens growth, but not on the rest. 

fit<-lm(mu.model~TestId*as.factor(concentration), data=parms.g, subset=parms.g$AddId=="Nicotine"|parms.g$AddId=="Control")
fit<-lm(mu.model~TestId*as.factor(concentration), data=parms.g, subset=parms.g$AddId=="Nicotine"|parms.g$AddId=="Control")

summary(fit)
drop1(fit, .~., test="F")
#

fit<-lm(mu.model~TestId+(concentration), data=parms.nop, 
        subset=parms.nop$AddId=="Catalpol"|parms.nop$AddId=="Control")
summary(fit)
drop1(fit, .~., test="F")
# yes, significant when all concentrations included
fit<-lm(mu.model~TestId+(concentration), data=parms.nop, 
        subset=parms.nop$AddId=="Catalpol"&parms.nop$concentration!=5|parms.nop$AddId=="Control")
summary(fit)
drop1(fit, .~., test="F")


fit<-lm(mu.model~TestId+(concentration), data=parms.nop, 
        subset=parms.nop$AddId=="Caffeine"|parms.nop$AddId=="Control")
summary(fit)
drop1(fit, .~., test="F")
fit<-lm(mu.model~TestId+(concentration), data=parms.nop, 
        subset=parms.nop$AddId=="Caffeine"&parms.nop$concentration!=5|parms.nop$AddId=="Control")
summary(fit)
drop1(fit, .~., test="F")


fit<-lm(mu.model~TestId+(concentration), data=parms.nop, 
        subset=parms.nop$AddId=="Nicotine"|parms.nop$AddId=="Control")
summary(fit)
drop1(fit, .~., test="F")
fit<-lm(mu.model~TestId+(concentration), data=parms.nop, 
        subset=parms.nop$AddId=="Nicotine"&parms.nop$concentration!=5|parms.nop$AddId=="Control")
summary(fit)


#########################

N<-levels(parms$AddId)
quartz()
par(mfrow=c(2,3))
for(i in 1:length(N)){
  m<-N[i]
  N.s<-parms[parms$AddId==m,]
  bargraph.CI(TestId, A.model, group=concentration,data=N.s, main=m)
}

#plot growth trials
# preprocessing code from (http://www.r-bloggers.com/analyzing-microbial-growth-with-r/)
####################
#install.packages("dplyr")
library(dplyr)
library(reshape2)
names(gr)[4]<-"Well"
annotated <- melt(gr, id=c("Species", "Nectar", "Concentration", "Well"), , 
                  variable.name="Time", value.name="OD595")
annotated$Time2<-as.numeric(substring(annotated$Time, 2, 5))
grouped <- group_by(annotated, Species, Nectar, Concentration)
stats <- summarise(grouped, N=length(OD595), Average=mean(OD595), StDev=sd(OD595))
conf_int95 <- function(data) {
    n <- length(data)
    error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
    return(error)
}
stats <- summarise(grouped, N=length(OD595), Average=mean(OD595),
                   CI95=conf_int95(OD595))
stats <- annotated %.%
          group_by(Species, Nectar, Concentration, Time2) %.%
          summarise(N=length(OD595), 
                    Average=mean(OD595, na.rm=T),
                    CI95=conf_int95(OD595)) %.%
          filter(!is.na(Concentration))


# plot
ggplot(data=stats, aes(x=Time2/60, y=Average, color=Species)) + 
       geom_line() + 
         facet_grid(Nectar ~ Concentration) +
       labs(x="Time (Hours)", y="Absorbance at 595 nm")

# plot and separate by Nectar and concentration
quartz()
ggplot(data=stats, aes(x=Time2/60, y=Average, color=Species)) + 
     geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=Species),
                   color=NA, alpha=0.3) + 
       geom_line() + 
       facet_grid(Nectar ~ Concentration) +
       labs(x="Time (Hours)", y="Absorbance at 595 nm")+
  ylim(0.1,0.4)

#pseudomonas only

stats.p<-stats[stats$Species=="P. fluorescens",]

quartz()
ggplot(data=stats.p, aes(x=Time2/60, y=Average, color=Species)) + 
       geom_line() + 
       facet_grid(Nectar ~ Concentration) +
       labs(x="Time (Hours)", y="Absorbance at 595 nm")+
  ylim(0.1,0.4)


ggplot(data=stats, aes(x=Time/3600, y=Average, color=Strain)) +
       geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=Strain),
                   color=NA, alpha=0.3) + 
       geom_line() +
       scale_y_log10() +
       facet_grid(Environment ~ .) +
       labs(x="Time (Hours)", y="Absorbance at 600 nm")
# add best-fit lines?

# plot only P.fluorescens
noP<-stats[stats$Species!="P. fluorescens",]

gr<-read.csv("rlv_microbe_gr1_2.csv", na.strings=c(".", "NA"), header=T)
gr.d<-data.frame(gr)

quartz()
ggplot(data=noP, aes(x=Time2/60, y=Average)) + 
       geom_line() + 
       facet_grid(Nectar ~ Concentration) +
       labs(x="Time (Hours)", y="Absorbance at 595 nm")

# plot from raw data
ggplot(data=grow.t)

n<-matrix(gr.d[,5:262])
ne<-as.data.frame(n[1,],times[1,])
plot(ne[,1]~ne[,2])


# figures for presentation (LSRF)
quartz()
parms.g$Microbe<-ifelse(parms.g$TestId=="A. astilbes"|parms.g$TestId=="Gluconobacter sp.", 1,2)
parms.g$Microbe[which(parms.g$TestId=="P. fluorescens")]<-3

quartz()
lineplot.CI(concentration, mu.model,group=TestId, data=parms.g, 
            subset=AddId=="Catalpol"&TestId!="P. fluorescens"|AddId=="Control"&TestId!="P. fluorescens", 
            legend=T, main="Catalpol", ylab="Max. growth rate (mu)", 
            xlab="Concentration", col=c(1,2,1,2,2), lwd=2, ylim=c(0,0.001))
lineplot.CI(concentration, mu.model,group=TestId, data=parms.g, 
            subset=AddId=="Catalpol"|AddId=="Control", 
            legend=T, main="Catalpol", ylab="Max. growth rate (mu)", 
            xlab="Concentration", col=c(1,2,1,2,3,2), lwd=2, ylim=c(0,0.001))

parms.nc<-parms[parms$TestId!="Control",]
quartz()
par(mfrow=c(1,2))
lineplot.CI(concentration, mu.model,group=TestId, data=parms.nc, 
            subset=AddId=="Nicotine"&TestId!="P. fluorescens"|AddId=="Control"&TestId!="P. fluorescens", 
            legend=T, main="Nicotine", ylab="Max. growth rate (mu)", 
            xlab="Concentration", col=c(1,2,1,2,2), lwd=2, ylim=c(0,0.001))
lineplot.CI(concentration, mu.model,group=TestId, data=parms.nc, 
            subset=AddId=="Nicotine"|AddId=="Control", 
            legend=T, main="Nicotine", ylab="Max. growth rate (mu)", 
            xlab="Concentration", col=c(1,2,1,2,3,2), lwd=2, ylim=c(0,0.001))
quartz()
par(mfrow=c(3,2))
lineplot.CI(concentration, mu.model,group=TestId, data=parms.nc, 
            subset=AddId=="Aucubin"|AddId=="Control", 
            legend=T, main="Aucubin", ylab="Max. growth rate (mu)", 
            xlab="Concentration", col=c(1,2,1,2,3,2), lwd=2, ylim=c(0,0.001))
lineplot.CI(concentration, mu.model,group=TestId, data=parms.nc, 
            subset=AddId=="Catalpol"|AddId=="Control", 
            legend=T, main="Catalpol", ylab="Max. growth rate (mu)", 
            xlab="Concentration", col=c(1,2,1,2,3,2), lwd=2, ylim=c(0,0.001))

lineplot.CI(concentration, mu.model,group=TestId, data=parms.nc, 
            subset=AddId=="Caffeine"|AddId=="Control", 
            legend=T, main="Caffeine", ylab="Max. growth rate (mu)", 
            xlab="Concentration", col=c(1,2,1,2,3,2), lwd=2, ylim=c(0,0.001))
lineplot.CI(concentration, mu.model,group=TestId, data=parms.nc, 
            subset=AddId=="Nicotine"|AddId=="Control", 
            legend=T, main="Nicotine", ylab="Max. growth rate (mu)", 
            xlab="Concentration", col=c(1,2,1,2,3,2), lwd=2, ylim=c(0,0.001))


lineplot.CI(concentration, mu.model,group=TestId, data=parms.g, 
            subset=AddId=="Ouabain"|AddId=="Control", 
            legend=T, main="Ouabain", ylab="Max. growth rate (mu)", 
            xlab="Concentration", col=c(1,2,1,2,3,2), lwd=2, ylim=c(0,0.001))
