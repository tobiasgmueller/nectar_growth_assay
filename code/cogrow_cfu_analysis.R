# this is the analysis for the co-growth experiments
# 5/16/21


# I grew 3 pairs of microbes across 6 different treatments for 3 days
# pair 1 = starmerella & zygosaccharomyces
# pair 2 = rosenbergiella and metschnikowia
# pair 4 = rosenbergiella and saccharomyces 

# the treatments were:
# 1 Control
# 2 Deltaline
# 3 Linalool
# 4 H202
# 5 Ethanol
# 6 tsb



#setwd("C:/Users/obiew/Desktop/vannette/growthassay_2021")


### packages ####

library(ggplot2)
library(dplyr)
library(gridExtra)
library(DescTools)
library(viridis)
library(reshape2)
library(ggpubr)
 
### read in data and clean ####

rm(list = ls()) 


p1 <- read.csv("input/cogrowthcfu_p1_notsb.csv")
p1<- na.omit(p1)
p1$pair <- as.factor(p1$pair)
p1$treatment <- as.factor(p1$treatment)


p2 <- read.csv("input/cogrowthcfu_p2_notsb.csv")
p2<- na.omit(p2)
p2$pair <- as.factor(p2$pair)
p2$treatment <- as.factor(p2$treatment)


p4 <- read.csv("input/cogrowthcfu_p4_notsb.csv")
p4<- na.omit(p4)
p4$pair <- as.factor(p4$pair)
p4$treatment <- as.factor(p4$treatment)



solo <- read.csv("input/cogrowthcfu_solo.csv")
solo<- na.omit(solo)
solo <- rename(solo, c("microbe" = "plate.ID"))
solo$microbe <- as.factor(solo$microbe)



# rosenbergiella is measured as a percent whereas yeasts are cfu
# i'll adjust rosenbergiella to be a percent of the max cfu 

p2$Rosenbergiella <- ((p2$Rosenbergiella/100)*max(p2$Metschnikowia))
p4$Rosenbergiella <- ((p4$Rosenbergiella/100)*max(p4$Saccharomyces))
solo$cfu[which(solo$microbe=="Rosenbergiella")] <- ((solo$cfu[which(solo$microbe=="Rosenbergiella")]/100)*max(solo$cfu))

# I should also adjust everything to the same dilution
# that is back calculate to cfu / ul solution

p1.adjust<-p1
p2.adjust<-p2
p4.adjust<-p4

p1.adjust$Zygosaccharomyces <- (p1$Zygosaccharomyces*20)
p1.adjust$Starmerella <- (p1$Starmerella*200)

p2.adjust$Rosenbergiella <- (p2$Rosenbergiella*2)
p2.adjust$Metschnikowia <- (p2$Metschnikowia*200)

p4.adjust$Rosenbergiella <- (p4$Rosenbergiella*2)
p4.adjust$Saccharomyces <- (p4$Saccharomyces*200)

solo$cfu.ul<- (solo$cfu*solo$dilution)
solo$treatment <- "solo"
solo$type <- "solo"


# raw cfu df
# dont use raw for analyses 
# use adjusted

melt1 <- p1[ -c(5) ]
melt1<-reshape2::melt(melt1, id=c("pair", "treatment"))
melt1 <- rename(melt1, c("microbe" = "variable"))
melt1 <- rename(melt1, c("cfu" = "value"))
melt1$type <- "co-inoculation"

s <- solo[solo$microbe == "Starmerella" | solo$microbe == "Zygosaccharomyces", ]
p1_cfu<- bind_rows(melt1, s)
p1_cfu<- p1_cfu[ -c(6:8) ]
p1_cfu <- droplevels(p1_cfu)

# and adjusted dilution df
melt1 <- p1.adjust[ -c(5) ]
melt1<-reshape2::melt(melt1, id=c("pair", "treatment"))
melt1 <- rename(melt1, c("microbe" = "variable"))
melt1 <- rename(melt1, c("cfu.ul" = "value"))
melt1$type <- "co-inoculation"

s <- solo[solo$microbe == "Starmerella" | solo$microbe == "Zygosaccharomyces", ]
p1_dil_adjust<- bind_rows(melt1, s)
p1_dil_adjust<- p1_dil_adjust[ -c(6:8) ]
p1_dil_adjust$treatment<-as.factor(p1_dil_adjust$treatment)
p1_dil_adjust <- droplevels(p1_dil_adjust)
p1_dil_adjust$treatment<- relevel(p1_dil_adjust$treatment, ref="solo")


### P2

melt2 <- p2[ -c(5) ]
melt2<-reshape2::melt(melt2, id=c("pair", "treatment"))
melt2 <- rename(melt2, c("microbe" = "variable"))
melt2 <- rename(melt2, c("cfu" = "value"))
melt2$type <- "co-inoculation"

s2 <- solo[solo$microbe == "Metschnikowia" | solo$microbe == "Rosenbergiella", ]
p2_cfu<- bind_rows(melt2, s2)
p2_cfu<- p2_cfu[ -c(6:8) ]
p2_cfu <- droplevels(p2_cfu)

# and adjusted dilution df
melt2 <- p2.adjust[ -c(5) ]
melt2<-reshape2::melt(melt2, id=c("pair", "treatment"))
melt2 <- rename(melt2, c("microbe" = "variable"))
melt2 <- rename(melt2, c("cfu.ul" = "value"))
melt2$type <- "co-inoculation"


s2 <- solo[solo$microbe == "Metschnikowia" | solo$microbe == "Rosenbergiella", ]
p2_dil_adjust<- bind_rows(melt2, s2)
p2_dil_adjust<- p2_dil_adjust[ -c(6:8) ]
p2_dil_adjust$treatment<-as.factor(p2_dil_adjust$treatment)
p2_dil_adjust <- droplevels(p2_dil_adjust)
p2_dil_adjust$treatment<- relevel(p2_dil_adjust$treatment, ref="solo")


### pair 4

melt4 <- p4[ -c(5) ]
melt4<-reshape2::melt(melt4, id=c("pair", "treatment"))
melt4 <- rename(melt4, c("microbe" = "variable"))
melt4 <- rename(melt4, c("cfu" = "value"))
melt4$type <- "co-inoculation"

s4 <- solo[solo$microbe == "Saccharomyces" | solo$microbe == "Rosenbergiella", ]
p4_cfu<- bind_rows(melt4, s4)
p4_cfu<- p4_cfu[ -c(6:8) ]
p4_cfu <- droplevels(p4_cfu)

# p4 adjusted
melt4 <- p4.adjust[ -c(5) ]
melt4<-reshape2::melt(melt4, id=c("pair", "treatment"))
melt4 <- rename(melt4, c("microbe" = "variable"))
melt4 <- rename(melt4, c("cfu.ul" = "value"))
melt4$type <- "co-inoculation"

s4 <- solo[solo$microbe == "Saccharomyces" | solo$microbe == "Rosenbergiella", ]
p4_dil_adjust<- bind_rows(melt4, s4)
p4_dil_adjust<- p4_dil_adjust[ -c(6:8) ]
p4_dil_adjust$treatment<-as.factor(p4_dil_adjust$treatment)
p4_dil_adjust <- droplevels(p4_dil_adjust)
p4_dil_adjust$treatment<- relevel(p4_dil_adjust$treatment, ref="solo")

#then lets also make dfs that dont have solo to compare to the control co growths
p1_dil_adjust2 <- p1_dil_adjust[p1_dil_adjust$treatment != "solo", ]
p1_dil_adjust2 <- droplevels(p1_dil_adjust2)

p2_dil_adjust2 <- p2_dil_adjust[p2_dil_adjust$treatment != "solo", ]
p2_dil_adjust2 <- droplevels(p2_dil_adjust2)

p4_dil_adjust2 <- p4_dil_adjust[p4_dil_adjust$treatment != "solo", ]
p4_dil_adjust2 <- droplevels(p4_dil_adjust2)

#------------------------- END OF DATA PREP



### Graphing ####


#for graphing we'll assign the same colors as our aggregate graphs 

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



pair1<- c("#882255","#44AA99")
pair2<- c("#117733","#DDCC77")
pair4<- c("#117733","#6699CC")
  

# graph pair 1
g1 <- ggplot(p1_cfu, aes(x=treatment, y=cfu, color=microbe)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_hline(yintercept=0)+
  labs(color = "Microbe pair")+
  ylab("Colony forming units")+
  theme_bw()+
  facet_grid(~type, scales = "free", space = "free")+
  theme(axis.title.x=element_blank(),
        legend.text = element_text(face = "italic"))


# dilution adjusted pair 1 log10

g2<- ggplot(p1_dil_adjust, aes(x=treatment, y=(log10(cfu.ul+1)),  color=microbe, fill=microbe)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_hline(yintercept=0)+
  labs(color = "Microbe pair", fill= "Microbe pair")+
  ylab("Log10 ( CFU/ul +1 )")+
  scale_fill_manual(values = pair1)+
  scale_color_manual(values = pair1)+
  theme_bw(base_size = 15)+
  scale_y_continuous(expand=expansion(mult = c(.1,.2)))+
  facet_grid(~type, scales = "free", space = "free")+
  theme(axis.title.x=element_blank(),
        legend.text = element_text(face = "italic"))


#grid.arrange(g1,g2)



# graph pair 2

g3 <- ggplot(p2_cfu, aes(x=treatment, y=cfu, color=microbe)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_hline(yintercept=0)+
  labs(color = "Microbe pair")+
  ylab("Colony forming units")+
  theme_bw()+
  facet_grid(~type, scales = "free", space = "free")+
  theme(axis.title.x=element_blank(),
        legend.text = element_text(face = "italic"))


g4 <- ggplot(p2_dil_adjust, aes(x=treatment, y=(log10(cfu.ul+1)), color=microbe, fill=microbe)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_hline(yintercept=0)+
  labs(color = "Microbe pair", fill= "Microbe pair")+
  ylab("Log10 ( CFU/ul +1 )")+
  theme_bw(base_size = 15)+
  scale_fill_manual(values = pair2)+
  scale_color_manual(values = pair2)+
  facet_grid(~type, scales = "free", space = "free")+
  scale_y_continuous(expand=expansion(mult = c(.1,.2)))+
  theme(axis.title.x=element_blank(),
        legend.text = element_text(face = "italic"))

#grid.arrange(g3,g4)


# graph pair 4

g5<- ggplot(p4_cfu, aes(x=treatment, y=cfu, color=microbe)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_hline(yintercept=0)+
  labs(color = "Microbe pair")+
  ylab("Colony forming units")+
  theme_bw()+
  facet_grid(~type, scales = "free", space = "free")+
  theme(axis.title.x=element_blank(),
        legend.text = element_text(face = "italic"))




g6 <- ggplot(p4_dil_adjust, aes(x=treatment, y=(log10(cfu.ul+1)), color=microbe, fill=microbe)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_hline(yintercept=0)+
  labs(color = "Microbe pair", fill= "Microbe pair")+
  ylab("Log10 ( CFU/ul +1 )")+
  theme_bw(base_size = 15)+
  scale_fill_manual(values = pair4)+
  scale_color_manual(values = pair4)+
  facet_grid(~type, scales = "free", space = "free")+
  scale_y_continuous(expand=expansion(mult = c(.1,.2)))+
  theme(axis.title.x=element_blank(),
        legend.text = element_text(face = "italic"))

#grid.arrange(g5,g6)





# graph just solo growth
g7<- ggplot(solo, aes(x=microbe, y=cfu)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  geom_hline(yintercept=0)+
  labs(color = "Microbe pair")+
  ylab("Colony forming units")+
  theme_bw()

#then adjust for dilutions
g8<- ggplot(solo, aes(x=microbe, y=cfu.ul)) +
  geom_boxplot(aes(), alpha=0.5) +
  geom_point(aes(), size=.5) +
  geom_hline(yintercept=0)+
  labs(color = "Microbe pair")+
  ylab("CFU adjusted")+
  theme_bw()





# just solo growth, non adjusted and adjusted
#grid.arrange(g7,g8)

# all graphs
#grid.arrange(g1, g2, g3, g4, g5, g6)

#non dilution corrected graphs
#ga1<- grid.arrange(g1,g3,g5)





# dilution adjusted graphs log transformed
# ga2<- grid.arrange(g2,g4,g6)
# above has legends on the sides which leads to uneven graph sizing but better layout


### FIGURE 4
ga2<- ggarrange(g2,g4,g6, 
                ncol=1, 
                nrow=3, 
                common.legend = FALSE, 
                legend="bottom",  
                labels=c("a","b","c"),
                label.x = 0.065,
                label.y = .985
                )
  
ga2
  ggsave(file="final_graphs/F4.svg", plot=ga2, width=180, height=270, units = "mm")
  ggsave(file="cogrowth_cfu_corrected.png", plot=ga2, width=8, height=8, units="in", dpi=300, )



  
  
### kruskal wallis analysis ####

# then a kruskal wallis to look for differences between treatments for microbes
# run seperately for each pair
kwtest.p1<- lapply(split(p1_dil_adjust, p1_dil_adjust$microbe), function(i){
  kruskal.test(cfu.ul ~ treatment, data = i)
})
kwtest.p1




kwtest.p1.cont<- lapply(split(p1_dil_adjust2, p1_dil_adjust2$microbe), function(i){
  kruskal.test(cfu.ul ~ treatment, data = i)
})
kwtest.p1.cont




kwtest.p2<- lapply(split(p2_dil_adjust, p2_dil_adjust$microbe), function(i){
  kruskal.test(cfu.ul ~ treatment, data = i)
})
kwtest.p2


kwtest.p4<- lapply(split(p4_dil_adjust, p4_dil_adjust$microbe), function(i){
  kruskal.test(cfu.ul ~ treatment, data = i)
})
kwtest.p4

#### dunnetts tests ####
# well they all come out significant. I think mostly becaues of the solo treatment where things grew

dunnets.cfu1<- lapply(split(p1_dil_adjust, p1_dil_adjust$microbe), function(i){
  DunnettTest(cfu.ul ~ treatment, data = i)
})

dunnets.cfu1.control<-lapply(split(p1_dil_adjust2, p1_dil_adjust2$microbe), function(i){
  DunnettTest(cfu.ul ~ treatment, data = i)
})


dunnets.cfu2<- lapply(split(p2_dil_adjust, p2_dil_adjust$microbe), function(i){
  DunnettTest(cfu.ul ~ treatment, data = i)
})

dunnets.cfu2

dunnets.cfu2.control<-lapply(split(p2_dil_adjust2, p2_dil_adjust2$microbe), function(i){
  DunnettTest(cfu.ul ~ treatment, data = i)
})


dunnets.cfu4<- lapply(split(p4_dil_adjust, p4_dil_adjust$microbe), function(i){
  DunnettTest(cfu.ul ~ treatment, data = i)
})

dunnets.cfu4


dunnets.cfu4.control<-lapply(split(p4_dil_adjust2, p4_dil_adjust2$microbe), function(i){
  DunnettTest(cfu.ul ~ treatment, data = i)
})


#### dunn test ####
# on second thought - I think much more interesting and related to our question
# is a pairwise comparing each treatment
# so dunn test with holm correction


dunn.cfu1<- lapply(split(p1_dil_adjust, p1_dil_adjust$microbe), function(i){
  DunnTest(cfu.ul ~ treatment, data = i, method="holm")
})
dunn.cfu1




dunn.cfu2<- lapply(split(p2_dil_adjust, p2_dil_adjust$microbe), function(i){
  DunnTest(cfu.ul ~ treatment, data = i, method="holm")
})

dunn.cfu2




dunn.cfu4<- lapply(split(p4_dil_adjust, p4_dil_adjust$microbe), function(i){
  DunnTest(cfu.ul ~ treatment, data = i, method="holm")
})


dunn.cfu4



