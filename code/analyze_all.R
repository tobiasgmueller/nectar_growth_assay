# analyze all
# this will aggregate the growth curves from the other r codes (bringing all plates together
# and analyze

# tobias mueller




# first load packages and all that good stuff
### packages ####

library(dplyr)
library(DescTools)
library(rstatix)
library(ggplot2)
library(reshape2)
library(lme4)
library(lmerTest)
library(bbmle)
library(emmeans) 
library(viridis)# for color pallettes
library(gridExtra) # ggplot editing
library(ggpubr) # ggplot editing
library(svglite) # to export graphs to SVGs
library(RColorBrewer)# for color pallettes
library(rcartocolor) # for color pallettes
library(FSA) # for dunn test
library(DHARMa) # for assumption testing
library(MASS)# for neg binomial
library(tidyverse)

rm(list = ls()) # cleans 

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


##### color pallette ####
#first lets make a good color scheme
safe_pal<-c("#DDCC77", #metsch
            "#CC6677", # acin
            "#117733", # rosen
            "#88CCEE", #panto
            "#332288", #aureo
            "#AA4499", # bacill
            "#44AA99", # starmerella
            "#999933",# rhodo
            "#888888", # pseudo
            "#661100", # pecto
            "#6699CC", # sacc
            "#882255") # zygo

### read in data and clean ####

parm_10xaa <- read.csv("output/highaa_summary.csv")
parm_10xaa$plate <- "10xaa"

parm_10xaa$treatment <- as.factor(parm_10xaa$treatment)
parm_10xaa$microbe <- as.factor(parm_10xaa$microbe)
parm_10xaa$plate <- as.factor(parm_10xaa$plate)
  

parm_2xaa <- read.csv("output/2xaa_summary.csv")
parm_2xaa$plate <- "2xaa"

parm_2xaa$treatment <- as.factor(parm_2xaa$treatment)
parm_2xaa$microbe <- as.factor(parm_2xaa$microbe)
parm_2xaa$plate <- as.factor(parm_2xaa$plate)



parm_2xsug <- read.csv("output/low_sugar_summary.csv")
parm_2xsug$plate <- "2xsug"

parm_2xsug$treatment <- as.factor(parm_2xsug$treatment)
parm_2xsug$microbe <- as.factor(parm_2xsug$microbe)
parm_2xsug$plate <- as.factor(parm_2xsug$plate)



parm_2mMh2o2 <- read.csv("output/2mM_h2o2_summary.csv")
parm_2mMh2o2$plate <- "2mM H2O2"

parm_2mMh2o2$treatment <- as.factor(parm_2mMh2o2$treatment)
parm_2mMh2o2$microbe <- as.factor(parm_2mMh2o2$microbe)
parm_2mMh2o2$plate <- as.factor(parm_2mMh2o2$plate)
  

parm_4mMh2o2 <- read.csv("output/4mM_h2o2_summary.csv")
parm_4mMh2o2$plate <- "4mM H2O2"

parm_4mMh2o2$treatment <- as.factor(parm_4mMh2o2$treatment)
parm_4mMh2o2$microbe <- as.factor(parm_4mMh2o2$microbe)
parm_4mMh2o2$plate <- as.factor(parm_4mMh2o2$plate)
  

parm_delta <- read.csv("output/high_deltaline_summary.csv")
parm_delta$plate <- "22 ug/ml deltaline"

parm_delta$treatment <- as.factor(parm_delta$treatment)
parm_delta$microbe <- as.factor(parm_delta$microbe)
parm_delta$plate <- as.factor(parm_delta$plate)


parm_delta31 <- read.csv("output/high_deltaline_summary31.csv")
parm_delta31$plate <- "22 ug/ml deltaline"

parm_delta31$treatment <- as.factor(parm_delta31$treatment)
parm_delta31$microbe <- as.factor(parm_delta31$microbe)
parm_delta31$plate <- as.factor(parm_delta31$plate)



parm_linalool <- read.csv("output/100ng_linalool_summary.csv")
parm_linalool$plate <- "100 ng/ml linalool"

parm_linalool$treatment <- as.factor(parm_linalool$treatment)
parm_linalool$microbe <- as.factor(parm_linalool$microbe)
parm_linalool$plate <- as.factor(parm_linalool$plate)

parm_ltp <- read.csv("output/150ng_LTP_summary.csv")
parm_ltp$plate <- "150 ug/ml LTP"

parm_ltp$treatment <- as.factor(parm_ltp$treatment)
parm_ltp$microbe <- as.factor(parm_ltp$microbe)
parm_ltp$plate <- as.factor(parm_ltp$plate)


parm_sug <- read.csv("output/2xsugar_summary.csv")
parm_sug$plate <- "30% sugar"

parm_sug$treatment <- as.factor(parm_sug$treatment)
parm_sug$microbe <- as.factor(parm_sug$microbe)
parm_sug$plate <- as.factor(parm_sug$plate)


parm_etoh <- read.csv("output/1perc_etoh_summary.csv")
parm_etoh$plate <- "1% EtOH"

parm_etoh$treatment <- as.factor(parm_etoh$treatment)
parm_etoh$microbe <- as.factor(parm_etoh$microbe)
parm_etoh$plate <- as.factor(parm_etoh$plate)



# then we want to aggregate these so we can get an overall control 
# make into 2 df, one for old nectar recipe, one for new recipe (10xaa)
# the old nectar recipe was worse (many things didnt grow in it) 
# but the data is included here incase

parm_all_old <- bind_rows(parm_10xaa, parm_2xaa, parm_2xsug)

parm_all <- bind_rows(parm_4mMh2o2, parm_2mMh2o2, parm_linalool,parm_delta, parm_ltp, parm_sug, parm_etoh)



#rename some columns left over from grofit package
parm_all <- rename(parm_all, c("date" = "concentration"))
parm_all <- rename(parm_all, c("strain" = "microbe"))

# reorder columns in df for readibility
parm_all <- parm_all[moveme(names(parm_all), "plate before date")]


# fix factors
parm_all$treatment <- as.factor(parm_all$treatment)
parm_all$strain <- as.factor(parm_all$strain)
parm_all$date <- as.factor(parm_all$date)  
parm_all$plate <- as.factor(parm_all$plate)  


# add microbe names
microbenames<- read.csv("input/microbe_names.csv")
parm_all <- merge(parm_all, microbenames, by = "strain", incomparables=NA)

# fix more factors
parm_all$microbe <- as.factor(parm_all$microbe)
parm_all$family <- as.factor(parm_all$family)
parm_all$order <- as.factor(parm_all$order)  
parm_all$class <- as.factor(parm_all$class)  
parm_all$kingdom <- as.factor(parm_all$kingdom)  


# and then also lets reorder the microbes from nectar specialists to non specialists
parm_all$microbe <- factor(parm_all$microbe, levels = c("Metschnikowia reukaufii", 
                                                        "Acinetobacter nectaris", 
                                                        "Rosenbergiella nectarea",
                                                        "Pantoea agglomerans" ,
                                                        "Aureobasidium pullulans" ,
                                                        "Bacillus subtilis",
                                                        "Starmerella bombi",
                                                        "Rhodotorula fujisanensis",
                                                        "Pseudomonas mandelii",
                                                        "Pectobacterium carotovorum",
                                                        "Saccharomyces cerevisiae",
                                                        "Zygosaccharomyces bailii"))



# now I'd like to adjust all values relative to control growth on each plate
# that is - how a microbes 2 controls grew relative to the mean controls of that microbe
# across all plates

# (mean of control / mean of control on plate) * each treatment value


# first calculate out the control means for each microbe on each plate 
  parm_all <- parm_all %>%
    group_by(plate, microbe) %>%
  mutate(control_mean_plate.A =  mean(A.model[which(treatment=="control")]))%>%
  ungroup()

  # and the overall control mean for each microbe
  parm_all <- parm_all %>%
    group_by(microbe) %>%
    mutate(control_mean.A =  mean(A.model[which(treatment=="control")]))%>%
    ungroup()
  

#then figure out the multiplier 
parm_all$multiplier.A <- parm_all$control_mean.A / parm_all$control_mean_plate.A
parm_all <- parm_all[moveme(names(parm_all), "multiplier.A before date")]

# okay lets now multiply all our Alpha values by their multiplier 
parm_all$adjusted.A <- parm_all$A.model*(parm_all$multiplier.A)
parm_all <- parm_all[moveme(names(parm_all), "adjusted.A before date")]


# then for ease of coding, i will replace A.model with the adjusted term
parm_all$A.model <- parm_all$adjusted.A


# now do the above all again for mu
# first calculate out the control means for each treatment
parm_all <- parm_all %>%
  group_by(plate, microbe) %>%
  mutate(control_mean_plate.mu =  mean(mu.model[which(treatment=="control")]))%>%
  ungroup()

# and the control mean for each microbe
parm_all <- parm_all %>%
  group_by(microbe) %>%
  mutate(control_mean.mu =  mean(mu.model[which(treatment=="control")]))%>%
  ungroup()


#then figure out the multiplier 
parm_all$multiplier.mu <- parm_all$control_mean.mu / parm_all$control_mean_plate.mu
parm_all <- parm_all[moveme(names(parm_all), "multiplier.mu before date")]

# okay lets now multiply all our Alpha values by their multiplier 
parm_all$adjusted.mu <- parm_all$mu.model*(parm_all$multiplier.mu)
parm_all <- parm_all[moveme(names(parm_all), "adjusted.mu before date")]


# then for ease of coding, i will replace mu.model with the adjusted term
parm_all$mu.model <- parm_all$adjusted.mu




#### Now to scale data so microbes can be compared to each other
# ill just adjust everything
# to be a percent of the control mean
# scaled from 0 (no growth), to 1 (equal to control / no change), to above 1 (increased growth)

parm_all <- parm_all %>%
  group_by(microbe) %>%
  mutate(scaled.A = (A.model/mean(A.model[which(treatment=="control")])))%>%
  ungroup()

parm_all <- parm_all[moveme(names(parm_all), "scaled.A before A.model")]



parm_all <- parm_all %>%
  group_by(microbe) %>%
  mutate(scaled.mu = (mu.model/mean(mu.model[which(treatment=="control")])))%>%
  ungroup()

parm_all <- parm_all[moveme(names(parm_all), "scaled.mu before mu.model")]



# and then im going to make a new column
# that puts controls for a microbe together, and treatments for each microbes
# this column will be called "type"

parm_all$type <- with(parm_all, ifelse (treatment=="control", 
                                        paste(parm_all$treatment), 
                                        paste(parm_all$plate)))

#then move to the front for easy seeing
parm_all <- parm_all[moveme(names(parm_all), "type before date")]
parm_all$type <- as.factor(parm_all$type)

#then relevel type so control is first. this is required for dunnetts test
parm_all$type<- relevel(parm_all$type, ref="control")

# we'll order treatment levels with control first
# and then most to least impactful
# so graphs are pretty

parm_all$type <- factor(parm_all$type, levels = c("control",
                                                  "4mM H2O2", 
                                                  "2mM H2O2", 
                                                  "100 ng/ml linalool",
                                                  "150 ug/ml LTP",
                                                  "30% sugar",
                                                  "1% EtOH",
                                                  "22 ug/ml deltaline"
                                                  ))

parm_all$plate <- factor(parm_all$plate, levels = c("control",
                                                  "4mM H2O2", 
                                                  "2mM H2O2", 
                                                  "100 ng/ml linalool",
                                                  "150 ug/ml LTP",
                                                  "30% sugar",
                                                  "1% EtOH",
                                                  "22 ug/ml deltaline"
))


# parm_all$type <- factor(parm_all$type, levels = c("control",
#                                                         "4mM H2O2", 
#                                                         "2mM H2O2", 
#                                                         "100 ng/ml linalool",
#                                                         "22 ug/ml deltaline" ,
#                                                         "150 ug/ml LTP",
#                                                         "30% sugar",
#                                                         "1% EtOH"))


#then also make a df of only treatment rows and only controls 
parm_treatonly <- subset(parm_all, treatment == "treatment")
parm_contonly <- subset(parm_all, treatment == "control")

# also make bacteria and yeast only df
parm_yeastonly <- subset(parm_all, kingdom == "yeast")
parm_bactonly <- subset(parm_all, kingdom == "Bacteria")




#_______________________________ END OF DATA PREP ________________________________

### Analysis begins____________________________________ ####




### ------------------ Krustak wallis ----------------------####

# this splits the df into a list of df based on microbe
# then runs a function (kruskal test) on each df in the list
# checks if within a microbe there are differences in parameter between treatments

kwtest.A<- lapply(split(parm_all, parm_all$microbe), function(i){
  kruskal.test(A.model ~ type, data = i)
})
kwtest.A


kwtest.mu<- lapply(split(parm_all, parm_all$microbe), function(i){
  kruskal.test(mu.model ~ type, data = i)
})
kwtest.mu





####  ------------------- Dunnettes test ----------------------------------#### 
#  run a dunnetts test by microbe 
# his computes pairwise differences between treatments and control for each microbe


dunnets.A<- lapply(split(parm_all, parm_all$microbe), function(i){
  DunnettTest(A.model ~ type, data = i)
})



dunnets.mu <- lapply(split(parm_all, parm_all$microbe), function(i){
  DunnettTest(mu.model ~ type, data = i)
})



#### ------------- negative binomial models ----------------------------------------

# to assess the overall ability for microbes to grow across treatments


nb.scaled.a<-glm.nb(data=parm_all, 
              scaled.A ~ type)

simoutput<-simulateResiduals(fittedModel = test)
plot(simoutput)



#test2<-glmer.nb(data=parm_all, 
#      scaled.A ~ type + (1|microbe))
# adding the random effect here doesnt change coefficients while increasing AIC



nb.scaled.mu<-glm.nb(data=parm_all, 
             scaled.mu ~ type)

test2<-glmer.nb(data=parm_all, 
                scaled.mu ~ type + (1|microbe))





##### yeast vs bacteria ####
#it could be fun to add something about kingdom to this analysis - look at yeasts vs bacteria
memk <- lmer(data=parm_treatonly, A.model ~ kingdom + (1|type)+ (1|microbe))
summary(memk)

AIC(memk)

simoutput<-simulateResiduals(fittedModel = memk)
plot(simoutput)

# than also look at scaled data
nb.k.a <- glmer.nb(data=parm_treatonly, scaled.A ~ kingdom + (1|type) + (1|microbe))
summary(nb.k.a)

# nb.k.a <- glm.nb(data=parm_treatonly, scaled.A ~ kingdom)
# summary(nb.k.a)

simoutput<-simulateResiduals(fittedModel = nb.k.a)
plot(simoutput)

# nb.k.mu <- glmer.nb(data=parm_treatonly, scaled.mu ~ kingdom + (1|type) + (1|microbe))
# summary(nb.k.mu)
# 
# nb.k.mu <- glm.nb(data=parm_treatonly, scaled.mu ~ kingdom )
# summary(nb.k.mu)





### levels comparison ####

# I want to run a comparison to see if there are trends
# between high / med / low nectar specialists


# first add a new column for specialization
parm_treatonly <- parm_treatonly %>% 
  mutate(rank = ifelse(microbe=="Metschnikowia reukaufii"|
                         microbe=="Acinetobacter nectaris"|
                         microbe=="Rosenbergiella nectarea"|
                         microbe=="Pantoea agglomerans", "high",
                       ifelse(microbe=="Aureobasidium pullulans"|
                                microbe=="Starmerella bombi"|
                                microbe=="Rhodotorula fujisanensis"|
                                microbe=="Bacillus subtilis"|
                                microbe=="Pseudomonas mandelii", "medium", "low")))

parm_treatonly$rank <- as.factor(parm_treatonly$rank)
parm_treatonly$rank <- factor(parm_treatonly$rank, levels= c("high","medium","low"))


# run kruskal wallis to test for differences in scaled Alpha and mu
kw.rank.a <- kruskal.test(data=parm_treatonly, scaled.A ~rank)
kw.rank.a # not significant difference

kw.rank.mu <- kruskal.test(data=parm_treatonly, scaled.mu ~rank)
kw.rank.mu # comes out significant

# followup dunntest on mu
DunnTest(scaled.mu~rank, data=parm_treatonly, method="holm")



# correlation test ####
cor.test(parm_treatonly$scaled.A,parm_treatonly$scaled.mu, method = "pearson",exact=FALSE)

# run a correlation for each microbe
corr.microbe<- lapply(split(parm_treatonly, parm_treatonly$microbe), function(i){
  cor.test(i$scaled.A, i$scaled.mu, method = "pearson",exact=FALSE)
})




### plots for final figures ####





# z score A faceted by treatment, color by microbe

g1 <- ggplot(parm_treatonly, aes(x=microbe, y=z.a, color=microbe, fill=microbe)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  geom_hline(yintercept=0)+ 
  ylab("Scaled mpact on max OD (sd)")+
  xlab("Treatment")+
  labs(fill ="Microbes", color="Microbes")+
  facet_wrap(~plate)+
  theme_bw()+
  scale_fill_manual(values = safe_pal)+
  scale_color_manual(values = safe_pal)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.text = element_text(face = "italic"))


ggsave(file="alphaz_by_treatment.svg", plot=g1, width=180, height=135, units = "mm")
ggsave(file="alphaz_by_treatment.png", plot=g1, width=8, height=6, units="in", dpi=300,)


##### Figure 1 ####
#pairwise comparisions for graphing
nb.A.pc <- emmeans(nb.scaled.a, pairwise~type, adjust="tukey", type="response") # type="respomse" is essential to make output not in log scale

treatmentMM <- data.frame(nb.A.pc$emmeans)
treatmentMM <- rename(.data = treatmentMM, plate=type, predicted=response,ymin=asymp.LCL, ymax=asymp.UCL)
treatmentMM <- treatmentMM[treatmentMM$plate!="control",]

treatmentMeans <- ggplot(parm_treatonly, aes(x=plate, y=scaled.A)) +
  #geom_boxplot(aes(), alpha=0.5) +
  geom_point(size=3,position= position_jitterdodge(.4), aes(color=microbe), alpha=.3) +
  geom_hline(yintercept=1)+ 
  geom_errorbar(data=treatmentMM,aes(x=plate, ymin=ymin, ymax=ymax),width=.3, inherit.aes = FALSE)+
  geom_point(data=treatmentMM, aes(y=predicted),color="black", fill="white", shape=21, size=2, stroke=1)+
  ylab("Scaled Impact on Max OD")+
  xlab("Treatment")+
  labs(color="Microbe", fill="Microbe")+
  theme_bw(base_size = 12)+
  scale_fill_manual(values = safe_pal)+
  scale_color_manual(values = safe_pal)+
  theme(axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.text = element_text(face = "italic"),
        axis.text.x = element_text(angle=90))

magic <- guides(colour = guide_legend(override.aes = list(alpha = 1))) 
treatmentMeans <- treatmentMeans + magic

treatmentMeans

ggsave(file="final_graphs/F1.svg", plot=treatmentMeans, width=180, height=135, units = "mm")
ggsave(file="final_graphs/F1.png", plot=treatmentMeans, width=8, height=6, units="in", dpi=300, )

dev.off()




##### figure 2####
g1_free<- ggplot(parm_treatonly, aes(x=plate, y=scaled.A, color=plate, fill=plate)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=1,
             position= position_jitterdodge(1),
             alpha=.5) +
  geom_hline(yintercept=1)+ 
  ylab("Scaled Impact on Max OD")+
  xlab("Treatment")+
  labs(fill ="Treatment", color="Treatment")+
  facet_wrap(~microbe, scales="free", 
             labeller = label_wrap_gen(width=10))+
  scale_fill_manual(values = safe_pal)+
  scale_color_manual(values = safe_pal)+
  scale_y_continuous(expand=expansion(mult = c(.1,.2)))+
  theme_bw(base_size = 12)+
  theme(strip.text = element_text(face = "italic"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.text = element_text(face = "italic"),
        legend.position = "bottom")



# scale_color_brewer(palette = "RdBu")+
#   scale_fill_brewer(palette = "RdBu")+

g1_free
ggsave(file="final_graphs/F2.svg", plot=g1_free, width=180, height=135, units = "mm")
ggsave(file="F2.png", plot=g1_free, width=8, height=6, units="in", dpi=300, )


##### figure 3 ####
g4<- ggplot(parm_treatonly, aes(x=plate, y=scaled.mu, color=plate, fill=plate)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  geom_hline(yintercept=1)+ 
  ylab("Scaled Impact on Growth Rate")+
  xlab("Treatment")+
  labs(color = "Treatment", fill="Treatment")+
  facet_wrap(~microbe, scale="free", 
             labeller = label_wrap_gen(width=10))+
  scale_y_continuous(expand=expansion(mult = c(.1,.2)))+
  theme_bw(base_size = 12)+
  scale_fill_manual(values = safe_pal)+
  scale_color_manual(values = safe_pal)+
  theme(strip.text = element_text(face = "italic"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.text = element_text(face = "italic"),
        legend.position = "bottom")
g4

ggsave(file="final_graphs/F3.svg", plot=g4,width=180, height=135, units = "mm")
ggsave(file="F3.png", plot=g4, width=8, height=6, units="in", dpi=300, )



#### figure S2 -correlation####
# facetted by microbe
p<- ggplot(parm_treatonly, 
       aes(x=scaled.mu, y=scaled.A))+
  geom_point(aes(color=plate), alpha=.5)+
  facet_wrap(~microbe, 
             labeller = label_wrap_gen(width=10))+
  geom_smooth(method="lm", 
              se=FALSE,
              fullrange=TRUE, 
              color="black",
              size=.7) +
  scale_fill_manual(values = safe_pal)+
  scale_color_manual(values = safe_pal)+
  theme(strip.text = element_text(face = "italic"))+ 
  xlab("Scaled Impact on Growth Rate")+
  ylab("Scaled Impact on Max OD")


# add text to each facet
dat_text <- data.frame(
  label = c("r = -0.50", "r = 0.93", "r = 0.17","r = 0.98","r = 0.13","r = 0.78","r = -0.44", "r = 0.36","r = 0.92","r = 0.98","r = 0.47","r = 0.63"
  ),
  microbe   = c("Metschnikowia reukaufii"  ,  "Acinetobacter nectaris" ,   
 "Rosenbergiella nectarea"  ,  "Pantoea agglomerans"   ,    
 "Aureobasidium pullulans"  ,  "Bacillus subtilis"     ,    
"Starmerella bombi"       ,   "Rhodotorula fujisanensis" , 
"Pseudomonas mandelii"    ,   "Pectobacterium carotovorum","Saccharomyces cerevisiae" ,  "Zygosaccharomyces bailii")  
)


corr_facet<- p + geom_text(
  data    = dat_text,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -0.1,
  vjust   = -9
)
corr_facet

ggsave(file="final_graphs/SF2.svg", plot=corr_facet,width=180, height=135, units = "mm")
ggsave(file="SF2.svg.png", plot=corr_facet, width=8, height=6, units="in", dpi=300, )



# all points together
corr_graph<- ggplot(parm_treatonly, 
       aes(x=scaled.mu, y=scaled.A))+
  geom_point(aes(color=plate), alpha=.5)+
  geom_smooth(method="lm", 
              se=FALSE,
              fullrange=TRUE, 
              color="black",
              size=.7) +
  scale_fill_manual(values = safe_pal)+
  scale_color_manual(values = safe_pal)+
  annotate(geom = "text", x = 2.2, y = 9, label = "r = 0.53 , p < .001")
corr_graph

ggsave(file="final_graphs/corr_graph.svg", plot=corr_graph,width=180, height=135, units = "mm")
ggsave(file="corr_graph.png", plot=corr_graph, width=8, height=6, units="in", dpi=300, )


# non scaled?
ggplot(parm_treatonly, 
       aes(x=mu.model, y=A.model))+
  geom_point(aes(color=plate), alpha=.5)+
  geom_smooth(method="lm", 
              se=FALSE,
              fullrange=TRUE, 
              color="black",
              size=.7) +
  scale_fill_manual(values = safe_pal)+
  scale_color_manual(values = safe_pal)+
  annotate(geom = "text", x = .002, y = 9, label = "Pearsons: 0.53 , P < .001")




# end of test zone



ggplot(parm_treatonly, aes(x=microbe, y=z.a, color=microbe, fill=microbe)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  geom_hline(yintercept=0)+ 
  ylab("Scaled impact on max OD (sd)")+
  xlab("Treatment")+
  labs(fill ="Microbes", color="Microbes")+
  facet_wrap(~plate, scales="free")+
  scale_y_continuous(expand=expansion(mult = c(.1,.2)))+
  theme_bw()+
  scale_fill_manual(values = bold)+
  scale_color_manual(values = bold)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.text = element_text(face = "italic"))




##### figure 4 ####
# must run correlation analysis code above to create rank factor
p1<-ggplot(parm_treatonly, aes(x=rank,y=scaled.A))+
  geom_boxplot()+
  geom_hline(yintercept=1)+
  theme_bw(base_size = 12)+ 
  ylab("Scaled Impact on Max OD")+
  xlab("Frequency of Isolation From Nectar")+
  scale_y_continuous(expand=expansion(mult = c(.1,.2)))

p2<-ggplot(parm_treatonly, aes(x=rank,y=scaled.mu))+
  geom_boxplot()+
  geom_hline(yintercept=1)+
  theme_bw(base_size = 12)+ 
  ylab("Scaled Impact on Growth Rate")+
  xlab("Frequency of Isolation From Nectar")+
  scale_y_continuous(expand=expansion(mult = c(.1,.2)))


freq_isolation<- ggarrange(p1,p2, labels=c("a","b"), hjust=-2)
freq_isolation
ggsave(file="final_graphs/F4.svg", plot=freq_isolation, width=180, height=100, units = "mm")
ggsave(file="F4.png", plot=freq_isolation, width=8, height=4, units="in", dpi=300, )






# quick graph to combine alpha and mu

microbe<- parm_treatonly$microbe
plate<- parm_treatonly$plate
a<- parm_treatonly$scaled.A
mu<- parm_treatonly$scaled.mu

dftest <- data.frame(microbe,plate,a,mu)
dftest <- melt(dftest, id.vars=c("microbe", "plate"))


levels(dftest$variable) <- c("Max OD","Growth rate")







##### figure S3 ####
# alpha non scaled faceted by microbe
g2 <- ggplot(parm_all, aes(x=type, y=A.model, color=type, fill=type)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  ylab("Maximum OD")+
  xlab("Treatment")+
  labs(fill ="Treatment", color="Treatment")+
  facet_wrap(~microbe, scales="free", 
             labeller = label_wrap_gen(width=10))+
  scale_fill_manual(values = safe_pal)+
  scale_color_manual(values = safe_pal)+
  scale_y_continuous(expand=expansion(mult = c(.1,.2)))+
  theme_bw(base_size = 12)+
  theme(strip.text = element_text(face = "italic"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.text = element_text(face = "italic"),
        legend.position = "bottom")


g2
ggsave(file="final_graphs/SF3.svg", plot=g2, width=180, height=135, units = "mm")
ggsave(file="SF3.png", plot=g2, width=8, height=6, units="in", dpi=300, )









##### figure S4 ####
# mu non scaled faceted by microbe
g5 <- ggplot(parm_all, aes(x=type, y=mu.model, color=type, fill=type)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  ylab("Growth Rate")+
  xlab("Treatment")+
  labs(color = "Treatment", fill="Treatment")+
  facet_wrap(~microbe, scale="free", 
             labeller = label_wrap_gen(width=10))+
  scale_y_continuous(expand=expansion(mult = c(.1,.2)))+
  theme_bw(base_size = 12)+
  scale_fill_manual(values = safe_pal)+
  scale_color_manual(values = safe_pal)+
  theme(strip.text = element_text(face = "italic", size=9),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.text = element_text(face = "italic"),
        legend.position = "bottom")
g5

ggsave(file="final_graphs/SF4.svg", plot=g5, width=180, height=135, units = "mm")
ggsave(file="SF4.png", plot=g5, width=8, height=6, units="in", dpi=300, )









# then just yeasts and just bacteria
##### figure S5 ####


# for some reason some fool (ahem toby..) made yeast lowercase
# so lets fix that
levels(parm_all$kingdom) <- list(Yeast  = "yeast", Bacteria = "Bacteria")

# okay i was real proud of my graph below with all the facets but its not 
# as relevant to the point so this is what were going with
sf5_kingdom<- parm_all %>%
  ggplot() +
  geom_boxplot(aes(x=kingdom, y=scaled.A, fill = kingdom),outlier.shape = NA, size=.7) +
  geom_jitter(aes(x=kingdom, y=scaled.A, color = microbe), size = 1.8, alpha = .5) +
  geom_hline(yintercept=1)+ 
  ylab("Scaled Impact on Max OD")+
  xlab("Kingdom")+
  labs(color = "Microbes", fill="Kingdom")+
  scale_y_continuous(expand=expansion(mult = c(.1,.2)))+
  theme_bw(base_size = 12)+
  scale_fill_manual(values = safe_pal)+
  scale_color_manual(values = safe_pal)+
  theme(axis.title.x=element_blank())
sf5_kingdom
ggsave(file="final_graphs/SF5.svg", plot=sf5_kingdom, width=180, height=135, units = "mm")




# and then also one thats facetted 
# itd be nice to have a grouped facet like margins() does in facet_grid

# so new function to achieve that
CreateAllFacet <- function(df, col){
  df$facet <- df[[col]]
  temp <- df
  temp$facet <- "all"
  merged <-rbind(temp, df)
  
  # ensure the facet value is a factor
  merged[[col]] <- as.factor(merged[[col]])
  
  return(merged)
}


df <- CreateAllFacet(parm_all, "type")

# the change to factor isnt working in the function so repeat that
df$facet<-as.factor(df$facet)

#then lets make levels same as other graphs
df$facet <- factor(df$facet, levels = c(
  "4mM H2O2", 
  "2mM H2O2", 
  "100 ng/ml linalool",
  "150 ug/ml LTP",
  "30% sugar",
  "1% EtOH",
  "22 ug/ml deltaline",
  "control",
  "all"
))

levels(df$kingdom) <- list(Yeast  = "yeast", Bacteria = "Bacteria")


sf5_all<-df %>%
  ggplot() +
  geom_boxplot(aes(x=kingdom, y=scaled.A, fill = kingdom),outlier.shape = NA, size=.7) +
  geom_jitter(aes(x=kingdom, y=scaled.A, color = microbe), size = 1.8, alpha = .5) +
  geom_hline(yintercept=1)+ 
  ylab("Scaled Impact on Max OD")+
  xlab("Kingdom")+
  labs(color = "Microbes", fill="Kingdom")+
  facet_wrap(~facet, scales="free")+
  scale_y_continuous(expand=expansion(mult = c(.1,.2)))+
  theme_bw(base_size = 12)+
  scale_fill_manual(values = safe_pal)+
  scale_color_manual(values = safe_pal)+
  theme(axis.title.x=element_blank())
sf5_all
ggsave(file="final_graphs/SF5_kingdom_all.svg", plot=sf5_all, width=200, height=135, units = "mm")

















#### GRAVEYARD #### 
# old code below


#### old sf5 ####
yeast_pallette <- c("#DDCC77", #metsch
                    "#332288", #aureo
                    "#44AA99", # starmerella
                    "#999933",# rhodo
                    "#6699CC", # sacc
                    "#882255") # zygo

bacteria_pallette <-c("#CC6677", # acin
                      "#117733", # rosen
                      "#88CCEE", #panto
                      "#AA4499", # bacill
                      "#888888", # pseudo
                      "#661100") # pecto







#adjust yeast and bacteria level order
parm_yeastonly$microbe <- factor(parm_yeastonly$microbe, levels = c("Metschnikowia reukaufii", 
                                                                    "Aureobasidium pullulans" ,
                                                                    "Starmerella bombi",
                                                                    "Rhodotorula fujisanensis",
                                                                    "Saccharomyces cerevisiae",
                                                                    "Zygosaccharomyces bailii"))

parm_bactonly$microbe <- factor(parm_bactonly$microbe, levels = c( 
  "Acinetobacter nectaris", 
  "Rosenbergiella nectarea",
  "Pantoea agglomerans" ,
  "Bacillus subtilis",
  "Pseudomonas mandelii",
  "Pectobacterium carotovorum"))



yeast <- ggplot(parm_yeastonly, aes(x=microbe, y=scaled.A, color=microbe, fill=microbe)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  geom_hline(yintercept=1)+ 
  ylab("Scaled Impact on Max OD")+
  xlab("Treatment")+
  ggtitle("Yeast")+
  labs(color = "Microbes", fill="Microbes")+
  facet_wrap(~type, scales="free")+
  scale_y_continuous(expand=expansion(mult = c(.1,.2)))+
  theme_bw(base_size = 12)+
  scale_fill_manual(values = yeast_pallette)+
  scale_color_manual(values = yeast_pallette)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.text = element_text(face = "italic"))


bacteria <- ggplot(parm_bactonly, aes(x=microbe, y=scaled.A, color=microbe, fill=microbe)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  geom_hline(yintercept=1)+ 
  ylab("Scaled Impact on Max OD")+
  xlab("Treatment")+
  ggtitle("Bacteria")+
  labs(color = "Microbes", fill="Microbes")+
  facet_wrap(~type, scales="free")+
  scale_y_continuous(expand=expansion(mult = c(.1,.2)))+
  theme_bw(base_size = 12)+
  scale_fill_manual(values = bacteria_pallette)+
  scale_color_manual(values = bacteria_pallette)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.text = element_text(face = "italic"))


g3 <- ggarrange(yeast, bacteria, nrow=1, legend="bottom")
g3
ggsave(file="final_graphs/SF5.svg", plot=g3, width=360, height=135, units = "mm")
ggsave(file="SF5.png", plot=g3, width=16, height=6, units="in", dpi=300, )












##### OLD -figure S2 ####
g_combo<- ggplot(dftest, aes(x=plate, y=value, fill=variable, color=variable))+
  geom_boxplot(aes(),  alpha=.5)+
  ylab("Scaled impact on growth rate and maximum OD") +
  xlab("Treatment")+
  geom_hline(yintercept=1, alpha=.6)+
  facet_grid(microbe~plate, scales = "free")+
  theme_bw(base_size = 14)+
  scale_y_continuous(expand=expansion(mult = c(.1,.1)))+
  theme(strip.text.y = element_text(angle=0),
        legend.title = element_blank()) 
g_combo



# new combo scatter plot
centroids <- aggregate(cbind(scaled.A,scaled.mu)~plate,parm_treatonly,mean)
f         <- function(z)sd(z)/sqrt(length(z)) # function to calculate std.err
se        <- aggregate(cbind(se.x=scaled.A,se.y=scaled.mu)~plate,parm_treatonly,f)
centroids <- merge(centroids,se, by="plate") 


g_combo2<- ggplot(parm_treatonly, aes(x=scaled.A, y=scaled.mu, color=factor(plate)))+
  geom_point(aes(),  alpha=.5)+ 
  geom_point(data=centroids, size=2)+
  geom_errorbar(data=centroids,aes(ymin=scaled.mu-se.y,ymax=scaled.mu+se.y),width=0.1)+
  geom_errorbarh(data=centroids,aes(xmin=scaled.A-se.x,xmax=scaled.A+se.x),height=0.1)+
  ylab("Scaled growth rate") +
  xlab("Scaled max OD")+
  geom_hline(yintercept=1, alpha=.6)+
  facet_wrap(~microbe, scales = "free")+
  theme_bw(base_size = 14)+
  scale_y_continuous(expand=expansion(mult = c(.1,.1)))+
  theme(strip.text.y = element_text(angle=0),
        legend.title = element_blank()) 
g_combo2

ggsave(file="final_graphs/SF2.svg", plot=g_combo, width=360, height=270, units = "mm")
ggsave(file="SF2.png", plot=g_combo, width=11, height=8, units="in", dpi=300, )


#### anova (ish) party ####

# and then an test  of microbe within treatment (only looking at treatment z scores)
# however the data is very non normal and the variance is quite different across groupings
# therefor 
# Kruskal wallis party!

#linalool
aov.lin <- aov(z.a ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "100 ng/ml linalool")  )
summary(aov.lin)
linpc <- emmeans(aov.lin, pairwise~microbe, adjust="tukey")
linpc


kw.lin<-kruskal.test(z.a ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "100 ng/ml linalool") )
kw.lin
FSA::dunnTest(z.a ~ microbe, 
              data=subset(parm_treatonly,parm_treatonly$plate == "100 ng/ml linalool"),
              method = "holm")

# 4mM h2o2
aov.h2o24 <- aov(z.a ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "4mM H2O2")  )
summary(aov.h2o24)
h2o24pc <- emmeans(aov.h2o24, pairwise~microbe, adjust="tukey")
h2o24pc


kw.h2o2<-kruskal.test(z.a ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "4mM H2O2")  )
kw.h2o2
FSA::dunnTest(z.a ~ microbe, 
              data=subset(parm_treatonly,parm_treatonly$plate == "4mM H2O2"),
              method = "holm")



# 2mM h2o2
aov.h2o22 <- aov(z.a ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "2mM H2O2")  )
summary(aov.h2o22)
h2o22pc <- emmeans(aov.h2o22, pairwise~microbe, adjust="tukey")
h2o22pc

#deltaline
aov.del <- aov(z.a ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "22 ug/ml deltaline")  )
summary(aov.del)
delpc <- emmeans(aov.del, pairwise~microbe, adjust="tukey")
delpc

#ltp
aov.ltp <- aov(z.a ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "150 ug/ml LTP")  )
summary(aov.ltp)
ltppc <- emmeans(aov.ltp, pairwise~microbe, adjust="tukey")
ltppc

# 2xsugar
aov.sug <- aov(z.a ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "30% sugar")  )
summary(aov.sug)
sugpc <- emmeans(aov.sug, pairwise~microbe, adjust="tukey")
sugpc


# 1% etoh
aov.etoh <- aov(z.a ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "1% EtOH")  )
summary(aov.etoh)
etohpc <- emmeans(aov.etoh, pairwise~microbe, adjust="tukey")
etohpc









# and then on non z scores


#linalool
aov.linA <- aov(A.model ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "100 ng/ml linalool")  )
summary(aov.linA)

# 4mM h2o2
aov.h2o24A <- aov(A.model ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "4mM H2O2")  )
summary(aov.h2o24A)

# 2mM h2o2
aov.h2o22A <- aov(A.model ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "2mM H2O2")  )
summary(aov.h2o22A)

#deltaline
aov.delA <- aov(A.model ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "22 ug/ml deltaline")  )
summary(aov.delA)

#ltp
aov.ltpA <- aov(A.model ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "150 ug/ml LTP")  )
summary(aov.ltpA)

# 2xsugar
aov.sugA <- aov(A.model ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "30% sugar")  )
summary(aov.sugA)

# 1% etoh
aov.etohA <- aov(A.model ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "1% EtOH")  )
summary(aov.etohA)

# for control
aov.contA<- aov(A.model ~ microbe, subset(parm_all,parm_all$treatment == "control")  )
summary(aov.contA)





#then look at mu

#control
aov.contmu<- aov(mu.model ~ microbe, subset(parm_all,parm_all$treatment == "control")  )
summary(aov.contmu)


aov.linmu <- aov(mu.model ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "100 ng/ml linalool")  )
summary(aov.linmu)

# 4mM h2o2
aov.h2o24mu <- aov(mu.model ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "4mM H2O2")  )
summary(aov.h2o24mu)

# 2mM h2o2
aov.h2o22mu <- aov(mu.model ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "2mM H2O2")  )
summary(aov.h2o22mu)

#deltaline
aov.delmu <- aov(mu.model ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "22 ug/ml deltaline")  )
summary(aov.delmu)

#ltp
aov.ltpmu <- aov(mu.model ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "150 ug/ml LTP")  )
summary(aov.ltpmu)

# 2xsugar
aov.sugmu <- aov(mu.model ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "30% sugar")  )
summary(aov.sugmu)

# 1% etoh
aov.etohmu <- aov(mu.model ~ microbe, subset(parm_treatonly,parm_treatonly$plate == "1% EtOH")  )
summary(aov.etohmu)





###### variance analysis ####
# and very much not equal variance. like wow. thats interesting




# then it could be interesting to compare the variance between each group

levene_test(z.a~rank,data=parm_treatonly)
levene_test(z.mu~rank,data=parm_treatonly)

library(onewaytests)
bf.test(z.a~rank,data=parm_treatonly)
bf.test(z.mu~rank,data=parm_treatonly)


# we can also just extract the residuals and run a pairwise test on those

#first calculate the median
treatonly_resid <- parm_treatonly %>% group_by(rank) %>%
  mutate(z.mu.med = median(z.mu)) %>%
  ungroup()

treatonly_resid <- treatonly_resid %>% group_by(rank) %>%
  mutate(z.a.med = median(z.a)) %>%
  ungroup()


#then the residuals
treatonly_resid$z.mu.res<-abs(treatonly_resid$z.mu-treatonly_resid$z.mu.med)
treatonly_resid$z.a.res<-abs(treatonly_resid$z.a-treatonly_resid$z.a.med)

ggplot(treatonly_resid, aes(x=rank,y=z.mu.res))+
  geom_boxplot()+
  theme_bw()+
  ylab("Scaled impact on max growth (sd)")+
  xlab("Frequency of isolation from nectar")

ggplot(treatonly_resid, aes(x=rank,y=z.a.res))+
  geom_boxplot()+
  theme_bw()+
  ylab("Scaled impact on max growth (sd)")+
  xlab("Frequency of isolation from nectar")


rankmuaov<-aov(data=treatonly_resid, z.mu.res~rank)
summary(rankmuaov)
TukeyHSD(rankmuaov)

rankmukw<-kruskal.test(data=treatonly_resid, z.mu.res~rank)
rankmukw
DunnTest(data=treatonly_resid, z.mu.res~rank, method="holm")


rankakw<-kruskal.test(data=treatonly_resid, z.a.res~rank)
rankakw
DunnTest(data=treatonly_resid, z.a.res~rank, method="holm")



##### correlation analysis ####

# okay first create a column "rank" ordered by nectar "specialization"
# really this is just abundance / frequency 
parm_treatonly$rank <- as.numeric(parm_treatonly$microbe)

# decided the above isnt great so im splitting microbes into 3 groups, 
# 1= "high", 2-"med", 3="low"


parm_treatonly <- parm_treatonly %>% 
  mutate(rank = ifelse(microbe=="Metschnikowia reukaufii"|
                         microbe=="Acinetobacter nectaris"|
                         microbe=="Rosenbergiella nectarea"|
                         microbe=="Pantoea agglomerans", 1,
                       ifelse(microbe=="Aureobasidium pullulans"|
                                   microbe=="Starmerella bombi"|
                                   microbe=="Rhodotorula fujisanensis"|
                                   microbe=="Bacillus subtilis"|
                                   microbe=="Pseudomonas mandelii", 2, 3)))



table(parm_treatonly$microbe,parm_treatonly$rank) # just quickly check that it worked

# then lets look for correlation across all plates
cor.test(parm_treatonly$rank,parm_treatonly$z.a, method = "spearman",exact=FALSE)


cor.test(parm_treatonly$rank,parm_treatonly$z.a, method = "kendall",exact=FALSE)



# and a matching graph
corg1<-ggplot(parm_treatonly, aes(rank, z.a, color=plate))+
  geom_point() +
  geom_smooth(method = "lm", formula=y~x)+
  ylab("Scaled impact on maximum growth")+
  scale_x_continuous(breaks = round(seq(min(parm_treatonly$rank), 
                                        max(parm_treatonly$rank), by =1),1),
                     labels = c("High","med","low"))+
  annotate("text", x = 1.5, y = 100,
           label = "Kendall correlation tau= -0.0003, p= .99")+
  theme(axis.text.x=element_text(angle = -70, hjust = 0),
        axis.title.x=element_blank())


# and for growth rate
cor.test(parm_treatonly$rank,parm_treatonly$z.mu, method = "spearman",exact=FALSE)
cor.test(parm_treatonly$rank,parm_treatonly$z.mu, method = "kendall",exact=FALSE)


corg2<- ggplot(parm_treatonly, aes(rank, z.mu, color=plate))+
  geom_point() +
  geom_smooth(method = "lm", formula=y~x)+
  ylab("Scaled impact on growth rate")+
  scale_x_continuous(breaks = round(seq(min(parm_treatonly$rank), 
                                        max(parm_treatonly$rank), by =1),1),
                     labels = c("High","med","low"))+
  annotate("text", x = 1.5, y = 100,
           label = "Kendall correlation tau= -0.13, p < .005")+
  theme(axis.text.x=element_text(angle = -70, hjust = 0),
        axis.title.x=element_blank())

ggarrange(corg1,corg2, common.legend = TRUE)

# and then run a correlation for each tretment individualy

# deltaline - mu is sig
delta <- subset(parm_treatonly, plate == "22 ug/ml deltaline")

cor.test(delta$rank,delta$z.a, method = "kendall",exact=FALSE)
cor.test(delta$rank,delta$z.mu, method = "kendall",exact=FALSE)


# 4mM h2o2
h2o24mm <- subset(parm_treatonly, plate == "4mM H2O2")

cor.test(h2o24mm$rank,h2o24mm$z.a, method = "kendall",exact=FALSE)
cor.test(h2o24mm$rank,h2o24mm$z.mu, method = "kendall",exact=FALSE)


# 2mM h2o2 - mu is sig
h2o22mm <- subset(parm_treatonly, plate == "2mM H2O2")

cor.test(h2o22mm$rank,h2o22mm$z.a, method = "kendall",exact=FALSE)
cor.test(h2o22mm$rank,h2o22mm$z.mu, method = "kendall",exact=FALSE)


# linalool
lina <- subset(parm_treatonly, plate == "100 ng/ml linalool")

cor.test(lina$rank,lina$z.a, method = "kendall",exact=FALSE)
cor.test(lina$rank,lina$z.mu, method = "kendall",exact=FALSE)

# LTP - mu is slightly sig but positive
ltp <- subset(parm_treatonly, plate == "150 ug/ml LTP")

cor.test(ltp$rank,ltp$z.a, method = "kendall",exact=FALSE)
cor.test(ltp$rank,ltp$z.mu, method = "kendall",exact=FALSE)


# 30% sugar -- mu is sig

sug <- subset(parm_treatonly, plate == "30% sugar")

cor.test(sug$rank,sug$z.a, method = "kendall",exact=FALSE)
cor.test(sug$rank,sug$z.mu, method = "kendall",exact=FALSE)


#1% ethanol
etoh <- subset(parm_treatonly, plate == "1% EtOH")

cor.test(etoh$rank,etoh$z.a, method = "kendall",exact=FALSE)
cor.test(etoh$rank,etoh$z.mu, method = "kendall",exact=FALSE)






#old graphs and old code graveyard






# z score A faceted by microbe, color by treatment
ggplot(parm_treatonly, aes(x=plate, y=z.a, color=plate)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  geom_hline(yintercept=0)+
  labs(color = "Treatment")+
  ylab("max OD z-score")+
  facet_wrap(~microbe)+
  theme_bw()+ 
  scale_color_viridis(discrete = TRUE)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        strip.text = element_text(face = "italic"))


#then alpha by absolute value (not z)
ggplot(parm_treatonly, aes(x=plate, y=A.model, color=plate)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  geom_hline(yintercept=0)+
  labs(color = "Treatment")+
  ylab("max OD")+
  facet_wrap(~microbe)+
  theme_bw()+
  scale_color_viridis(discrete = TRUE)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        strip.text = element_text(face = "italic"))


# and then mu
ggplot(parm_treatonly, aes(x=plate, y=z.mu, color=plate)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  geom_hline(yintercept=0)+
  labs(color = "Treatment")+
  ylab("growth rate z-score")+
  facet_wrap(~microbe)+
  theme_bw()+
  scale_color_viridis(discrete = TRUE)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        strip.text = element_text(face = "italic"))



# graphing controls to look for plate effects

ggplot(parm_contonly, aes(x=plate, y=A.model, color=microbe)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  theme(axis.text.x=element_text(angle = -70, hjust = 0)) +
  labs(title = "Alpha controls") +
  facet_wrap(~ microbe, scales = "free")



ggplot(parm_contonly, aes(x=plate, y=z.a)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  theme(axis.text.x=element_text(angle = -70, hjust = 0)) +
  labs(title = "Alpha z controls")



ggplot(parm_contonly, aes(x=plate, y=z.a, color=microbe)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  theme(axis.text.x=element_text(angle = -70, hjust = 0)) +
  labs(title = "alpha controls")







####------------------- Graphing!  ----------------------####

# look at how each treatment impacts by microbe (z score)
ggplot(parm_all, aes(x=type, y=z.a, color=treatment)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  theme(axis.text.x=element_text(angle = -70, hjust = 0)) +
  labs(title = "alpha z") +
  facet_grid(~ microbe, scale="free", space="free_x")


#then look mircobial growth by treatment (z score)
ggplot(parm_treatonly, aes(x=microbe, y=z.a, color=microbe)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  geom_hline(yintercept=0)+ 
  ylab("Carry capacity z-score")+
  xlab("treatment")+
  labs(color = "Microbes")+
  facet_wrap(~ plate)+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.text = element_text(face = "italic"))



# then lets also just do z scores by treatment for all microbes
ggplot(parm_treatonly, aes(x=plate, y=z.a, color=plate)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  geom_hline(yintercept=0)+
  labs(color = "Treatment")+
  ylab("Carry capacity z-score")+
  facet_wrap(~ microbe)+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        strip.text = element_text(face = "italic"))


# it should also be noted that the degree to which microbes can grow
# at all in nectar is vastly different
# below is plot of control OD values
ggplot(parm_contonly, aes(x=microbe, y=A.model, color=microbe)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  theme(axis.text.x=element_text(angle = 70, hjust = 1)) +
  labs(title = "alpha of control")

ggplot(parm_contonly, aes(x=microbe, y=mu.model, color=treatment)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  theme(axis.text.x=element_text(angle = 70, hjust = 1)) +
  labs(title = "mu z")


#  also quickly double check that control z score center acound zero
ggplot(parm_contonly, aes(x=type, y=z.a, color=treatment)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  theme(axis.text.x=element_text(angle = -70, hjust = 0)) +
  labs(title = "alpha z") +
  facet_grid(~ microbe, scale="free", space="free_x")

# and double check that the sd of these is 1
parm_contonly %>%
  group_by(microbe) %>%
  summarise(sd = sd(z.a))
# okay they all have a sd of 1, thats good

parm_contonly %>%
  group_by(microbe) %>%
  summarise(mean = mean(z.a))
# hmm but their means arnt all 0 - is that not how its supposed to be?



test <- by(parm_contonly$z.a,parm_contonly$microbe,mean)


#check for normality of z scores
ggplot(parm_all, aes(x=z.a))+
  geom_density()
# dont seem normal to me but I should probably run some more tests there



# then just some quick summaries of the data
z.a.summary<- parm_all %>%
  group_by(microbe, type) %>%
  summarise(
    count = n(),
    mean = mean(z.a, na.rm = TRUE),
    sd = sd(z.a, na.rm = TRUE),
    median = median(z.a, na.rm = TRUE),
    IQR = IQR(z.a, na.rm = TRUE)
  )%>%
  ungroup()

a.summary<- parm_all %>%
  group_by(microbe, type) %>%
  summarise(
    count = n(),
    mean = mean(A.model, na.rm = TRUE),
    sd = sd(A.model, na.rm = TRUE),
    median = median(A.model, na.rm = TRUE),
    IQR = IQR(A.model, na.rm = TRUE)
  ) %>%
  ungroup()







