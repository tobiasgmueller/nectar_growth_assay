# ------ Title ---
# Data analysis with Growthcurver

# 18 feb, 2020
#Tobias Mueller

#---  load packages -----
library(growthcurver)
library(lubridate)
library(magrittr)
library(reshape2)
library(ggplot2)
library(readr)
library(dplyr)
library(broom)
library(emmeans)
library(grofit)


# read in data

setwd("C:/Users/obiew/Desktop/vannette/96wellplate_growthassay")

# --- for plate #1 ec52 M. reukaufii - feburary 11th
# d <- read.csv("input/growth_feb11_600.csv")
# labels <- read.csv("input/well_labels.csv")

# --- for plate #2 ec102 A. pallulans (weird growth) - feburary 14th
#d <- read.csv("input/growth_feb14_600.csv")
#labels <- read.csv("input/well_labels.csv")

# --- restructed plate layout below this

# --- for rerun of ec52 - march 19th - plate is a little weird (i messed uo) so has own labels
d <- read.csv("input/ec52_19march2020.csv")
d_nocont <- read.csv("input/ec52_19march2020_nocontrol.csv")
labels <- read.csv("input/labels_16mar2020.csv")
labels_no_control <- read.csv("input/labels_16mar2020_nocontrol.csv")


#for rerun of ec102 A. pallulans - june 11th
# d <- read.csv("input/nectarchem_ec102_11jun2020_tm.csv")
# d_nocont <- d %>% select(-A7, -B7, -C7, -D1, -D2,
#                           -D3, -D4, -D5, -D7, -D8,
#                           -D9, -D10, -D11, -D12,
#                           -E6,-F6, -G6, -H6)
# labels <- read.csv("input/labels.csv")
# labels_no_control <- read.csv("input/labels_nocontrol.csv")



#for bb362 run on june 17th
# d <- read.csv("input/nectar_chem_bb362_17jun2020_tm.csv")
# d_nocont <- d %>% select(-A7, -B7, -C7, -D1, -D2,
#                          -D3, -D4, -D5, -D7, -D8,
#                          -D9, -D10, -D11, -D12,
#                          -E6,-F6, -G6, -H6)
# labels <- read.csv("input/labels.csv")
# labels_no_control <- read.csv("input/labels_nocontrol.csv")

#for q1f2 run on june 23th
# d <- read.csv("input/nectar_chem_q1f2_23jun2020_tm.csv")
# d_nocont <- d %>% select(-A7, -B7, -C7, -D1, -D2,
#                          -D3, -D4, -D5, -D7, -D8,
#                          -D9, -D10, -D11, -D12,
#                          -E6,-F6, -G6, -H6)
# labels <- read.csv("input/labels.csv")
# labels_no_control <- read.csv("input/labels_nocontrol.csv")


#for ec84 run on june 29th
# d <- read.csv("input/nectar_chem_ec84_29jun2020_tm.csv")
# d_nocont <- d %>% select(-A7, -B7, -C7, -D1, -D2,
#                          -D3, -D4, -D5, -D7, -D8,
#                          -D9, -D10, -D11, -D12,
#                          -E6,-F6, -G6, -H6)
# labels <- read.csv("input/labels.csv")
# labels_no_control <- read.csv("input/labels_nocontrol.csv")







#incase you need to check class type
#sapply(d, class)

# below is only importatant if treatment doesnt read as class factor
#labels$treatment<- as.factor(labels$treatment)


d <- rename(d, c("time" = "Time"))
d_nocont <- rename(d_nocont, c("time" = "Time"))
# if you get an error here that should be fine and just means the columns are correctly labled as "time"


#simplify time points to single hours - I dont think the package can deal with time - atleast I havnt gotten it to...
d$time <- hms(d$time)
d$time <- hour (d$time)

d_nocont$time <- hms(d_nocont$time)
d_nocont$time <- hour (d_nocont$time)




####### EXAMPLE with a single growth curve # -----
## 
#following: https://cran.r-project.org/web/packages/growthcurver/vignettes/Growthcurver-vignette.html

#for simplicity
d2 <- d_nocont





#---- removing wells with no growth ------


#this creates a new df that lists the percent change of each well then removes ones that are less than a cutoff (5% here  i.e  (x > .05))
# change the value of v (below) to change what the cutoff point is
v<-.05

#the first line is writing a function that calculates the percent change between the initial value and the value n rows below it
pct <- function(x) {x/lag(x, n=288)-1}
#then applying that to our raw data and removing things under out threshold
temp <- d[-1] %>%
  mutate_all(funs(pct)) %>%
  na.omit()%>% 
  select_if(function(x) any (x > v))

#then because im not too sure how pipelines im switchgin to base r to just merge it back with the original df
temp <- bind_rows(temp,d)
#then drop columns that didnt have the growth required
temp <- temp %>%
  select_if(~ !any(is.na(.)))
#then drop first row (percent growth)
growthwells <- temp[-1,]
#now add back time - in this case in minutes -- though technically the first read is 13min 38 seconds not 15 minutes -- I may fix this later
growthwells$time <- seq.int(nrow(growthwells))*15
#then put time at the from of growthwells
growthwells<- growthwells[,c(51,1:50)]
#im sure it could be done cleaner but it works







# ----- full plate pdf curves -----


# begin by Analyzeing entire plate
gc_out_noc <- SummarizeGrowthByPlate(d2, bg_correct ="min", plot_fit = TRUE, plot_file="fitplot_nocont.pdf")
#gc_out2 <- SummarizeGrowthByPlate(dat2, bg_correct ="min", plot_fit = TRUE, plot_file="fitplot.pdf")
#gc_out_noc <- merge(gc_out_noc, labels_no_control, by.x="sample", by.y="well", incomparables=NA)
head(gc_out_noc)


gc_out <- SummarizeGrowthByPlate(d, bg_correct ="min", plot_fit = TRUE, plot_file="fitplot.pdf")
#gc_out2 <- SummarizeGrowthByPlate(dat2, bg_correct ="min", plot_fit = TRUE, plot_file="fitplot.pdf")
#gc_out <- merge(gc_out, labels, by.x="sample", by.y="well", incomparables=NA)
head(gc_out)


gc_out_cut <- SummarizeGrowthByPlate(growthwells, bg_correct ="min", plot_fit = TRUE, plot_file="fitplot_cut.pdf")


# ---- DATA CHOP ----
#from here on out, only wells that met the growth criteria are being used




# ---- ggplot prep----- 
#this will crash if you run with wells that dont have curves modled to them (wells with no growth)
#the minimum growth criteria should fix this everytime now though

###### choose how to cut data ######

#df <- growthwells
df <- d2



#write all model data to a df
models.all <- lapply(df[2:ncol(df)], function(x) SummarizeGrowth(df$time, x))

#create adjusted OD by well over time
df.predicted.plate <- data.frame(time = df$time)
for (i in names(df[2:ncol(df)])) 
{df.predicted.plate[[i]] <- (models.all[[i]]$data$N)}

#create predicted OD by well over time
df.predicted.plate2 <- data.frame(time = df$time)
for (i in names(df[2:ncol(df)])) 
{df.predicted.plate2[[i]] <- predict(models.all[[i]]$model)}




#grill that cheese
melt1 <- melt(df, id.vars = "time", variable.name = "sample")
              #we dont want to use "od" here from original df as it hasnt been adjusted like the predicted od has
           
melt2 <- melt(df.predicted.plate, id.vars = "time", variable.name = "sample", value.name = "od")
melt3 <- melt(df.predicted.plate2, id.vars = "time", variable.name = "sample", value.name = "pred.od")
df.final <- cbind(melt1, od=melt2[,3], pred.od=melt3[,3])

#df.final2 <- merge(melt2, melt3)

rm(melt1)
rm(melt2)

#merge labels with data  - (without controls)
labels_no_control$treatment <- as.factor(labels_no_control$treatment)
labels_no_control$well <- as.factor(labels_no_control$well)
df.final2 <- merge(df.final, labels_no_control, by.x="sample", by.y="well", incomparables=NA)





# ---------------- graph growth curves -------------------

#first individual growth curves, overlayed by treatment (only growth wells)------ 
ggplot(df.final2, aes(x=time, y=od)) +
  geom_point(aes(), alpha=0.5) +
  geom_point(aes(y=pred.od), color="red", size=.5) +
  facet_wrap(df.final2$treatment) +
  theme_bw()


#graph all wells growth curves, in rows by treatment (only growth wells)
ggplot(df.final2, aes(x=time, y=od)) +
  geom_point(aes(), alpha=0.5) +
  geom_point(aes(y=pred.od), color="red", size=.5) +
  facet_grid(vars(df.final2$treatment), vars(replicate)) +
  theme_bw()









# merge labels with gc_out
#head(labels)
#head(gc_out)
gc_out_lab <- merge(gc_out, labels, by.x="sample", by.y="well", incomparables=NA)
gc_out_lab_noc <- merge(gc_out_noc, labels, by.x="sample", by.y="well", incomparables=NA)
gc_out_lab_cut <- merge(gc_out_cut, labels, by.x="sample", by.y="well", incomparables=NA)
# check that merge worked
#dim(labels)
#dim(gc_out)
#dim(gc_out_lab)
#head(gc_out_lab)


#---------- select here which dataframes you want to use for graphing and analysis---------------
#write the one you want to df

#cut well that didnt meet growth minimum
df <- gc_out_lab_cut

#no wells cut
#df <- gc_out_lab

#cut out wells with no microbes added
#df <- gc_out_lab_noc


# make long dataframe for graphing
dfmelt <-melt(df[,c("sample","treatment", "r","k", "t_mid", "auc_l", "auc_e", "t_gen")])



#now to reorder treatments
df$treatment <- factor(df$treatment, levels = c(
  "control no microbe",
  "control with microbe", 
  "high sugar", 
  "low sugar", 
  "high aa", 
  "low aa",
  "high H2O2",
  "low H2O2",
  "high LTP",
  "low LTP",
  "high linalool",
  "low linalool"
))

dfmelt$treatment <- factor(dfmelt$treatment, levels = c(
  "control no microbe",
  "control with microbe", 
  "high sugar", 
  "low sugar", 
  "high aa", 
  "low aa",
  "high H2O2",
  "low H2O2",
  "high LTP",
  "low LTP",
  "high linalool",
  "low linalool"
))

#get rid of control no microbe as its making y axis crazy
temp<-filter(dfmelt, treatment != c("NA"))
melt_growth2<-filter(temp, treatment != c("control no microbe"))

#do the same with gc out lab
temp<-filter(df, treatment != c("NA"))
gc_out_lab2<-filter(temp, treatment != c(
  "control_no_microbe"))


# --------------------------------- Plots! -------------------------

# r = max growth rate
# k = max carrying capacity 
# n0 = original population number
# t_gen = doubling time
# t_mid = time to inflection point (when pop = 1/2 k)
# auc_l = area under curve (smoothed logorithmic curve)
# auc_e = area under curve (experimental data not smoothed)

#graph all meta data by treatment !
ggplot(melt_growth2, aes(x=treatment, y=value, color = variable))+
  geom_point()+
  geom_boxplot()+
  facet_grid(melt_growth2$variable, scales="free")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))




#graph less of everything 
temp <- filter(melt_growth2, variable %in% c("r","k","t_mid"))

ggplot(temp, 
  aes(x=treatment, y=value, color = variable))+
  geom_point()+
  geom_boxplot()+
  facet_grid(temp$variable, scales="free")+
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))




# --------- specific graphs (zoomed in basically) -------------


#graphing r
ggplot(filter(gc_out_lab, treatment != c("control no microbe")),
       aes(x=treatment, y=r, color= treatment))+
  geom_point()+
  geom_boxplot()+
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))

qqnorm(gc_out_lab2$r)
qqline(gc_out_lab2$r)


#graphing k 
ggplot(gc_out_lab2, aes(x=treatment, y=k, color= treatment))+
  geom_point()+
  geom_boxplot()+
  theme_bw() + 
#  ylim(0,5) +
  theme(axis.text.x = element_text(angle = 90))
   

qqnorm(gc_out_lab2$k)
qqline(gc_out_lab2$k)


#graph inflection point time
ggplot(gc_out_lab2, aes(x=treatment, y=t_mid, color= treatment))+
  geom_point()+
  geom_boxplot()+
  theme_bw()

#graph gen  time
ggplot(gc_out_lab2, aes(x=treatment, y=t_gen, color= treatment))+
  geom_point()+
  geom_boxplot()+
  theme_bw()







# #------------ run big anovas across treatments--------
# for r 
 aovr <- aov(r ~ treatment, data=gc_out_lab)
TukeyHSD(aovr)
 

# for k 
 aovk <- aov(k ~ treatment, data=gc_out_lab)
 TukeyHSD(aovk)

 #for auc_l

 
 aovaucl <- aov(auc_l ~ treatment, data=gc_out_lab)
 TukeyHSD(aovaucl)
 
 
 
 
# ----- looking at sugar only (k) -----------
ggplot(data = filter(gc_out_lab, treatment %in% c(
  "control with microbe", "high sugar", "low sugar")), 
       aes(x=treatment, y=k))+
  geom_point()+
  geom_boxplot()+
  theme_bw()

aovksugar <- aov(k ~ treatment, data = filter(gc_out_lab, treatment %in% c(
  "control with microbe", "high sugar", "low sugar")))

TukeyHSD(aovksugar)









# ------- looking at AA only (k) ---------
ggplot(filter(gc_out_lab, treatment %in% c(
  "control with microbe", "high aa", "low aa")),
  aes(x=treatment, y=k))+
  geom_point()+
  geom_boxplot()+
  theme_bw()


aovkaa <- aov(k ~ treatment, data = filter(gc_out_lab, treatment %in% c(
  "control with microbe", "high aa", "low aa")))

TukeyHSD(aovkaa)


# ------- looking at AA only (r) ---------
ggplot(filter(gc_out_lab, treatment %in% c(
  "control with microbe", "high aa", "low aa")),
  aes(x=treatment, y=r))+
  geom_point()+
  geom_boxplot()+
  theme_bw()


aovkaa <- aov(r ~ treatment, data = filter(gc_out_lab, treatment %in% c(
  "control with microbe", "high aa", "low aa")))

TukeyHSD(aovkaa)





#------looking at LTP only ------
ggplot(filter(gc_out_lab, treatment %in% c(
  "control with microbe", "high LTP", "low LTP")),
  aes(x=treatment, y=k))+
  geom_point()+
  geom_boxplot()+
  theme_bw()


aovkLTP <- aov(k ~ treatment, data = filter(gc_out_lab, treatment %in% c(
  "control with microbe", "high LTP", "low LTP")))

TukeyHSD(aovkLTP)






# ----- looking at H2O2 only -----

ggplot(filter(gc_out_lab, treatment %in% c(
  "control with microbe", "high H2O2", "low H2O2")),
  aes(x=treatment, y=k))+
  geom_point()+
  geom_boxplot()+
  theme_bw()


aovkh2o2 <- aov(k ~ treatment, data = filter(gc_out_lab, treatment %in% c(
  "control with microbe", "high H2O2", "low H2O2")))

TukeyHSD(aovkh2o2)






# --- looking at linalool only ----
ggplot(filter(gc_out_lab, treatment %in% c(
  "control with microbe", "high linalool", "low linalool")),
  aes(x=treatment, y=k))+
  geom_point()+
  geom_boxplot()+
  theme_bw()

aovklina <- aov(k ~ treatment, data = filter(gc_out_lab, treatment %in% c(
  "control with microbe", "high linalool", "low linalool")))

TukeyHSD(aovklina)