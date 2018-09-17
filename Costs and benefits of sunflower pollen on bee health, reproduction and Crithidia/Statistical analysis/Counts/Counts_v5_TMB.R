#Analysis of Crithidia counts in microcolonies fed different pollens
#Jessica Leslie honors thesis 2015
#updated 2016.02.29 to remove "BeeID" random effect 
#(model was over-specified)

#script updated 2016.02.16 to use updated file "updated survival"
#mistakes where there were 2 rows for a single bee are now removed
#some wing measures have been added

#EPY: Updated 2017.06.03 to use glmmTMB model fit

rm(list = ls())
#load packages
#Libraries:
##Load packages 
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

#set working directory
rm(list=ls())
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/jess stats/Counts/R_analyses")
Data<-read.csv("C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/Pollen consumption trials/updated survival.csv")#rename some columns
View(Data)

#check dataset
str(Data)

#eliminate a priori microcolonies 12, 16 (uninfected treatment found infected)
# and 50, 54 (died within first week, did not collect data)
Data<- subset(Data, !MC==12 & !MC==16 & !MC==50 & !MC==54)
summary(Data$MC)

#rename some columns
Data$Microcolony<-Data$MC
Data$Microcol<-factor(Data$Microcol) #converts to factor
Data$Microcol #76 levels, good
Data$Colony<-Data$innoc..Date #inoculation date identifies which colony was used
Data$Pollen<-Data$pollen.treat
Data$daysinfect<-Data$days.infect
Data$totcrith<-Data$crithidia.count
#Data$totcrith<-Data$crithidia.count #converts count to cells/uL (option)
Data$wingmm<-(13/20)*Data$marginal.cell #converts wing measurement from ocular units to mm

Data
str(Data)

#now calculate size dimorphism for each colony
library(plyr)
#first make data frame with max and min size for each bee in the colony
#example:
#ddply(dt,~group,summarise,mean=mean(age),sd=sd(age))
DM<-ddply(Data,~Microcol, summarise,bigbee=max(wingmm),smallbee=min(wingmm))
str(DM) #good, data frame
#calculate colony dimorphism: ratio of largest to smallest bee minus 1
DM$dimorph<-(DM$bigbee/DM$smallbee) -1
DM
#merge this back to original data:
Data.DM<-join(x=Data, y=DM, by=c("Microcol"), type="left")
Data.DM
View(Data.DM)
length(Data$Microcol)#380
length(Data.DM$Microcol)#380, good
#use only infected bees in counts analysis:
Infected<-subset(Data.DM, !infect.treat=="U")
Infected$infect.treat #some problems here, "I" recognized as 2 levels
Infected<-droplevels(Infected)
str(Infected$infect.treat)
Infected$infect.treat
library(plyr)
Infected$infect.treat<-revalue(Infected$infect.treat, c("I "="I"))#merge the 2 levels
str(Infected$infect.treat) #good, down to 1 level
#eliminate bees without counts
Counted<-subset(Infected, !totcrith=="NA")


str(Counted$Pollen) #2 levels
Data<-Counted #renames
Data$BeeID<-Data$Bee.ID #renaming
Data$colony<-Data$Colony #renaming
Data<-droplevels(Data) #remove empty factor levels
Data<-dplyr::rename(Data, Count = totcrith)
############     MODELING        #######
#$Count as response variable
#$pollen as fixed effect, 
#Covariates of $daysinfect, $dimorph
#random effects: (blocking term accounting for non-independence of bees within microcolony within colony)
# + 1( colony|Microcol|BeeID )
#A model of nested random e???ects (block within site)
#would be 1|site/block; 
#a model of crossed random e???ects (block and year) would be:
#(1|block)+(1|year).

#negative binomial
#Start with full nesting
Winged<-subset(Data, !wingmm=="NA" ) #remove some NA's
Dimorphed<-subset(Winged, !dimorph=="NA") #remove some NA's

nbinom0<-glmmTMB(totcrith ~ Pollen  +  wingmm + dimorph + (1 |Colony/ Microcol), 
                  family = "nbinom2",
                  data=Dimorphed)  
summary(nbinom0) #remove dimorphism
nbinom01<-glmmTMB(totcrith ~ Pollen  + wingmm + (1 |Colony/ Microcol), 
                  family = "nbinom2",
                  data=Winged)  
summary (nbinom01) #let's remove wing size
nbinom001<-glmmTMB(totcrith ~ Pollen  +  (1 |Colony/ Microcol), 
                   family = "nbinom2",
                   data=Data) #great, no error!
summary(nbinom001) #big effect of pollen
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   3.6897     0.4030   9.155  < 2e-16 ***
#   PollenS      -2.2311     0.4988  -4.473  7.7e-06 ***

#Test significance: 
drop1(nbinom001, test = "Chisq")
# Df    AIC    LRT  Pr(>Chi)    
# <none>    597.86                     
# Pollen  1 617.34 21.478 3.579e-06 ***
##Ooh significant

CountFinal<-nbinom001

################### Plots #############
Lsm<-lsmeans(CountFinal, ~Pollen)
Pairs<-contrast(Lsm, "pairwise") 
Pairs
# contrast estimate        SE df z.ratio p.value
# B - S    2.231108 0.4987569 NA   4.473  <.0001

##Sunflowers distinguish themselves

##Get letters to plot
Lets<-cld(Lsm, Letters = letters)
Lets
Lets$.group<-gsub(" ", "", Lets$.group) ##remove spaces

################## PLOT ON SCALE OF RESPONSE ##########
#Let's exponentiate original mean+SE and mean-SE
Lsmeans.df<-as.data.frame(summary(lsmeans(CountFinal, ~Pollen)))
Lsmeans.df$Count<-exp(Lsmeans.df$lsmean)
Lsmeans.df$btlo<-exp(Lsmeans.df$lsmean - Lsmeans.df$SE) #lower bound
Lsmeans.df$bthi<-exp(Lsmeans.df$lsmean + Lsmeans.df$SE) #upper bound

Lsmeans.df

#Merge with letters
Lets2<- dplyr::select(Lets, Pollen, .group)
Lets2

Lsmeans.df.merged<- join(Lsmeans.df, Lets2, by = "Pollen")

Lsmeans.df.merged
# Pollen   lsmean        SE df asymp.LCL asymp.UCL     Count      btlo      bthi .group
# 1      B 3.689676 0.4030217 NA 2.8997682  4.479584 40.031882 26.753211 59.901280      b
# 2      S 1.458568 0.5066856 NA 0.4654823  2.451653  4.299797  2.590581  7.136721      a
Lsmeans.df<-Lsmeans.df.merged #overwrites old

#x-axis:
pollen.axis<-c("Buck", "Sun")
#y-axis:
ylabel<- expression(paste("Parasite load (cells * 0.02 ", mu, L^-1, "))", sep="") ) #check for encoding

#ready to plot?
p<-ggplot(Lsmeans.df, aes(x = Pollen, y = Count)) + 
  #Option to sprinkle raw data... 
  geom_point(data = Data, aes( x = Pollen, y = Count), position = "jitter", 
             color = "blue", shape = 25, size = 0.5, alpha = 0.4)+
  geom_pointrange(aes(ymin = btlo, ymax = bthi), size = 1, color = "black") +
  scale_x_discrete(labels=pollen.axis) +
  ylab(ylabel) + #y label
  xlab("Pollen")
p
#Add text
Panel.1d <- p + geom_text(aes(y = bthi + 0.4 * max(bthi), label = .group),
                          size = 8)
Panel.1d
# ggsave("Microcolony.counts.v3.TMB.pdf", height = 5, width =5)

