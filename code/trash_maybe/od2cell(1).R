# this will be a script to convert OD to cell counts

# 5/25/21
# Tobias Mueller




# first prep the environment
### packages ####
setwd("C:/Users/obiew/Desktop/vannette/growthassay_2021")

library(dplyr)

rm(list = ls())

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


# then read in the data

labels <- read.csv("input/ODtocell_labels.csv")
d <- read.csv("input/cell_to_OD_test_dilutions.csv")


# transpose df
d <- t(d)
d <- as.data.frame(d)

# get the row labels as column 1
setDT(d, keep.rownames = TRUE)[]
# make time the headers

# redataframe it
d <- as.data.frame(d)


#fix df formate
names(d) <- d[1,]
d <- d[-1,]


# now compute the average OD between the two reads
colnames(d)[1] <- "wells"
colnames(d)[2] <- "read1"
colnames(d)[3] <- "read2"

d$wells<-as.factor(d$wells)
d$read1<-as.numeric(d$read1)
d$read2<-as.numeric(d$read2)

d$od <- rowMeans(d[,c('read1', 'read2')], na.rm=TRUE)

# now add dilution data to OD values
d$dilution <- ifelse(grepl("A", d$wells), ".25", 
                     ifelse(grepl("B", d$wells), ".25", 
                            ifelse(grepl("C", d$wells), ".125", 
                                   ifelse(grepl("D", d$wells), ".125",
                                          ifelse(grepl("E", d$wells), ".0625", 
                                                 ifelse(grepl("F", d$wells), ".0625", 
                                                        ifelse(grepl("G", d$wells), ".3333", 
                                                               ifelse(grepl("H", d$wells), ".3333", "no"))))))))


d$dilution<-as.numeric(d$dilution)

# now merge those OD values with the labels csv to get OD for different cellcountes
curve <- merge(d, labels, by="wells")

# now get rid of extra columns we dont need and correct names
curve = subset(curve, select = -c(dilution.x,row,microbe) )
curve <- rename(curve, c("dilution" = "dilution.y"))





# now a quick test plot


ggplot(curve, aes(x=cells_per_ul, y=od)) +
  geom_point(aes(), size=.5) +
  geom_smooth(method=lm)+
  ylab("OD")+
  xlab("cells per ul")+
  labs(color = "Microbes")+
  facet_wrap(~ Microbe, scales="free")+
  theme_bw()

# ooof


ggplot(curve, aes(x=dilution, y=od)) +
  geom_point(aes(), size=.5) +
  geom_smooth(method=lm)+
  ylab("OD")+
  xlab("dilution")+
  labs(color = "Microbes")+
  facet_wrap(~ Microbe, scales="free")+
  theme_bw()



# now what if having the two curves together is throwing things off
# so lets split thw two dilution sets

curve_e <- curve %>%
  dplyr::filter(grepl('A|C|E|G', wells))

curve_f <- curve %>%
  dplyr::filter(grepl('B|D|F|H', wells))




p1<-ggplot(curve_e, aes(x=dilution, y=od)) +
  geom_point(aes(), size=.5) +
  geom_smooth(method=lm)+
  ylab("OD")+
  xlab("dilution")+
  labs(color = "Microbes")+
  facet_wrap(~ Microbe, scales="free")+
  theme_bw()

p2<-ggplot(curve_f, aes(x=dilution, y=od)) +
  geom_point(aes(), size=.5) +
  geom_smooth(method=lm)+
  ylab("OD")+
  xlab("dilution")+
  labs(color = "Microbes")+
  facet_wrap(~ Microbe, scales="free")+
  theme_bw()


grid.arrange(p1,p2)
