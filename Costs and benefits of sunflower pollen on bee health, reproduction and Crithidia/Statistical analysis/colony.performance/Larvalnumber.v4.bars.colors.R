####2016.03#########
#analysis of sunflower microcolony performance

#EPY: Updated 2017.06.03 to use glmmTMB model fit
#edited 2017.07.01: bar chart, resort significance letters

#####Larval production################
####Source data: use summary created by LS Adler
#EPY added back colonies from which new adults emerged -- 
#colonies 3, 11, 23, 25, 33
#coded as "1" for pupal production

#rm(list=ls()) #clear memory

#Load packages


setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/jess stats/colony.performance")

#Libraries:
##Load packages 
library(MASS)
library(glmmTMB)
library(car)
library(bbmle) #Aic tables
library(lsmeans) #post-hoc tests with finite size correction
library(multcomp) #post-hoc tests
library(plyr)
library(dplyr)
source("glmmTMB.to.lsmeans.R") #helper function for TMB to talk to lsmeans
library(ggplot2)
library(cowplot) #multi panel plots


Perf0<- read.csv('C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/jess stats/colony.performance/performance_sum_SAS_CALLOWS.ADDED.csv',
                 na.strings=".")
Bees<-read.csv("C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/Pollen consumption trials/updated survival.csv")
#use "bees" to get a column for dimorphism
Data<-Bees
Data$wingmm<-(13/20)*Data$marginal.cell #converts wing measurement from ocular units to mm
#first make data frame with max and min size for each bee in the colony
#example:
Data$Microcolony<-Data$MC #renaming
library(plyr)
DM<-ddply(Data,~Microcolony, summarise,bigbee=max(wingmm),smallbee=min(wingmm))
str(DM) #good, data frame
#calculate colony dimorphism: ratio of largest to smallest bee minus 1
DM$dimorph<-(DM$bigbee/DM$smallbee) -1
DM
#merge this back to original data:
Perf0$Microcolony<-Perf0$microcolony
Data.DM<-join(x=Perf0, y=DM, by=c("Microcolony"), type="left")
Data.DM$Microcolony
Perf0<-Data.DM

#familiarize with the Performance data
str(Perf0)
#exclude microcolonies 12, 16 (uninfected treatment, found infected), 
#also exclude 50, 54 (died within first week)
Perf02<-subset(Perf0, !Microcolony=="50" & !Microcolony=="54" & !Microcolony=="12" & !Microcolony=="16")
Perf02$Microcolony<-as.factor(Perf02$Microcolony) #changes from numeric to factor
str(Perf02$Microcolony) #check that it worked
str(Perf02)
str(Perf02$sourcecol) #good , 4 levels
library(dplyr)
Perf02$sourcecol
#View(Perf02)
#rename(Perf02, Colony=sourcecol) #not working
Perf02$Colony<-Perf02$sourcecol
str(Perf02) #duplicated, renamed the column
Perf03<-Perf02
str(Perf03$infect) #good , 2 levels
Perf03$Infection<-Perf03$infect
Perf03$Infection<-relevel(Perf03$Infection, ref="U")
#reference level now "U"
str(Perf03) #renamed the column

Perf03$Pollen<-Perf03$pollen

#scout the variable
Nothing<-subset(Perf03, larvwt=NA)
NoL <- Perf03[is.na(Perf03$larvwt),]
NoL$larvwt
length(NoL$larvwt)
xtabs(~Pollen, Perf03)
xtabs(~Pollen, NoL) #21 buckwheat, 14 sunflower
xtabs(~Infection, NoL) #17 uninfected, 18 infected
xtabs(~Colony, NoL) #JML-3 and JML-4 fewer larvae, >50% no larvae


Makers<-subset(Perf03, !larvwt=="NA")
Makers$larvwt
length(Makers$larvwt) #41 colonies w/ larvae

#LSA recommends analyzing larval count and larval mass separately
hist(Perf03$larvnum)
LCOUNTS<-ddply(Perf03,~Infection*Pollen, 
          summarise,numlarvae=mean(larvnum),
          meanmass=mean(larvwt, na.rm=TRUE))
LCOUNTS

#plenty of zeroes
#check distribution
qqp(Perf03$larvnum, "norm") #poor fit for low and high values

#try lognormal
# lnorm means lognormal
qqp(Perf03$larvnum, "lnorm") #excellent
#log-normal would be a biologically appropriate distribution
#negative binomial: first fit distribution parameters, then plot
nbinom <- fitdistr(Perf03$larvnum, "Negative Binomial")
qqp(Perf03$larvnum, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
nbinom <- fitdistr(1 + Perf03$larvnum, "Negative Binomial")
qqp(1 + Perf03$larvnum, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
#looks great

#try poisson
poisson <- fitdistr(Perf03$larvnum, "Poisson")
qqp(Perf03$larvnum, "pois", poisson$estimate)
#too many zeroes?, poor fit to full distrib...
#but maybe OK when accounting for treatment
Buck<-subset(Perf03, Pollen=="B")
Sun<-subset(Perf03, Pollen=="S")

poisson <- fitdistr(Buck$larvnum, "Poisson")
qqp(Buck$larvnum, "pois", poisson$estimate)
#not very good
poisson <- fitdistr(Sun$larvnum, "Poisson")
qqp(Sun$larvnum, "pois", poisson$estimate)
#again not very good

#try gamma
gamma <- fitdistr(Perf03$larvnum, "gamma")
gamma <- fitdistr(1+Perf03$larvnum, "gamma")
qqp(1+ Perf03$larvnum, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])
#OK, but this gamma is for continuous distribution

#let's compare Poisson and negative binomial models
library(glmmTMB)
#first try poisson model:
fitPoiss<- glmmTMB(larvnum ~ Infection*Pollen  + (1|Colony), 
                    family = "poisson",
                    data=Perf03)
summary(fitPoiss)
# AIC      BIC   logLik deviance df.resid 
# 467.0    478.7   -228.5    457.0       71 
#Deviance >> Residuals.... --> Overdispersed

fitNB<- glmmTMB(larvnum ~ Infection*Pollen  + (1|Colony), 
                    family = "nbinom2",
                    data=Perf03)
summary(fitNB) #fixed effects same
library(bbmle)
plot(residuals(fitPoiss))
plot(residuals(fitNB)) #not much difference
AICtab(fitPoiss, fitNB)
# dAIC  df
# fitNB      0.0 6 
# fitPoiss 132.3 5 
#big edge to nbinom

summary(fitNB)
drop1(fitNB, test = "Chisq")
# Df    AIC     LRT Pr(>Chi)
# <none>              334.67                 
# Infection:Pollen  1 332.82 0.14738   0.7011

#drop interaction
NB1<-update(fitNB, ~. - Infection:Pollen)
drop1(NB1, test = "Chisq")
# Df    AIC    LRT Pr(>Chi)   
# <none>       332.82                   
# Infection  1 331.57 0.7442 0.388320   
# Pollen     1 339.19 8.3728 0.003809 **

###CHOSE THIS NB1 AS FINAL MODEL########

################### Plots #############
Lsm<-lsmeans(NB1, ~Pollen)
Pairs<-contrast(Lsm, "pairwise") 
Pairs
# contrast  estimate        SE df z.ratio p.value
# B - S    -1.094325 0.3629958 NA  -3.015  0.0026

##Sunflowers distinguish themselves

##Get letters to plot
Lets<-cld(Lsm, Letters = letters, sort = FALSE) #alphabetical left to right
Lets
Lets$.group<-gsub(" ", "", Lets$.group) ##remove spaces

################## PLOT ON SCALE OF RESPONSE ##########
#Let's exponentiate original mean+SE and mean-SE
Lsmeans.df<-as.data.frame(summary(lsmeans(NB1, ~Pollen)))
Lsmeans.df$Count<-exp(Lsmeans.df$lsmean)
Lsmeans.df$btlo<-exp(Lsmeans.df$lsmean - Lsmeans.df$SE) #lower bound
Lsmeans.df$bthi<-exp(Lsmeans.df$lsmean + Lsmeans.df$SE) #upper bound

Lsmeans.df

#Merge with letters
Lets2<- dplyr::select(Lets, Pollen, .group)
Lets2

Lsmeans.df.merged<- join(Lsmeans.df, Lets2, by = "Pollen")

Lsmeans.df.merged
# Pollen    lsmean        SE df  asymp.LCL asymp.UCL    Count      btlo     bthi .group
# 1      B 0.1888019 0.6013676 NA -0.9898568  1.367461 1.207802 0.6619498 2.203770      a
# 2      S 1.2831265 0.5791604 NA  0.1479930  2.418260 3.607902 2.0217553 6.438444      b
#3x as many larvae

Lsmeans.df<-Lsmeans.df.merged #overwrites old

#x-axis:
pollen.axis<-c("Buck", "Sun")
#y-axis:
ylabel<-"Number of Larvae" 

#ready to plot?
pd<-position_dodge(width=0.9)
p<-ggplot(Lsmeans.df, aes(x = Pollen, y = Count)) + 
  geom_bar(aes(fill = Pollen), color = "black", 
           stat = "identity", position = pd)+
  geom_errorbar(aes(ymin = btlo, ymax = bthi), 
                size = 1, width=0.4, color = "black", position = pd) +
  scale_x_discrete(labels=pollen.axis) +
  ylab(ylabel) + #y label
  xlab("Pollen")+
  theme(legend.position = "none")+
  theme(axis.line = element_line(size = 2),
        text = element_text (face = "bold"))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.2 * max(Lsmeans.df$bthi)))+
  scale_fill_manual(values = Colors)
p
#Add text
Lnum<-p + geom_text(aes(y = bthi + 0.1 * max(bthi), label = .group),
                        size = 8)
Lnum

##Export
# setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/jess stats/colony.performance")
# ggsave("LARVALNUMBER.v5.reletter.pdf", height = 5, width =5)

xtabs( ~Pollen, Makers)
# Pollen
# B  S 
# 15 26 
