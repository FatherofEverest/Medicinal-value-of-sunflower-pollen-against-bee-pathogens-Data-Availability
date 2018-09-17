#load packages
library(nlme)
library(lme4)
library(car)
library(plyr)
library(lsmeans)
library(ggplot2)

rm(list=ls())
setwd("C:/Users/Evan/Dropbox/Jess Sum15/jess stats/Consumption")
Data<-read.csv(
  "C:/Users/Evan/Dropbox/Jess Sum15/Pollen consumption trials/updated corrected consumption values only.csv",
  na.strings = ".")
View(Data)
str(Data)
Data$Colony<-Data$Start.Day
Data$Time<-Data$Day
Data$Infection<-Data$Infect.
#change reference level
library(plyr)
Data$Infection <- relevel(Data$Infection, ref = "U")

Data$Microcolony<-as.factor(Data$Micro...)
#eliminate colonies 12 and 16, these were uninfected treatments that had infection.
#also eliminate colonies 50 and 54-- died within first week
Data0<-subset(Data, !Microcolony==12 & !Microcolony==16 & !Microcolony==50 & !Microcolony==54)
Data$MicrocolonySize<-Data$X..of.bees
Data$NRate<-Data$Nectar.per.bee.per.hour
Data$PRate<-Data$Corrected.pollen.per.bee.per.hour
str(Data)

###add size dimorphism (optional)
#I prefer not to do have this in model because there are so many values missing
SizeData<-read.csv("C:/Users/Evan/Dropbox/Jess Sum15/Pollen consumption trials/updated survival.csv")#rename some columns
#eliminate a priori microcolonies 12, 16 (uninfected treatment found infected)
# and 50, 54 (died within first week, did not collect data)
SizeData<- subset(SizeData, !MC==12 & !MC==16 & !MC==50 & !MC==54)
summary(SizeData$MC)

#rename some columns
SizeData$Microcolony<-SizeData$MC
SizeData$Microcol<-factor(SizeData$Microcol) #converts to factor
SizeData$Microcol #76 levels, good
SizeData$Colony<-SizeData$innoc..Date #inoculation date identifies which colony was used
SizeData$Pollen<-SizeData$pollen.treat
SizeData$daysinfect<-SizeData$days.infect
SizeData$totcrith<-SizeData$crithidia.count
#SizeData$totcrith<-SizeData$crithidia.count #converts count to cells/uL (option)
SizeData$wingmm<-(13/20)*SizeData$marginal.cell #converts wing measurement from ocular units to mm

SizeData
str(SizeData)

#now calculate size dimorphism for each colony
library(plyr)
#first make SizeData frame with max and min size for each bee in the colony
#example:
#ddply(dt,~group,summarise,mean=mean(age),sd=sd(age))
DM<-ddply(SizeData,~Microcol, summarise,bigbee=max(wingmm),smallbee=min(wingmm))
DM$Microcolony<-DM$Microcol
str(DM) #good, data frame
#calculate colony dimorphism: ratio of largest to smallest bee minus 1
DM$dimorph<-(DM$bigbee/DM$smallbee) -1
DM
#merge this back to original data:
Data.DM<-join(x=Data, y=DM, by=c("Microcolony"), type="left")

Data<-Data.DM
Data<-subset (Data, !PRate=="NA") #NA's caused errors in lme model

#####data inspection (plotting)
#see some superb examples of plotting by Ben Bolker
#look for his "Lab 1 mixed models"
#http://ms.mcmaster.ca/~bolker/classes/uqam/mixedlab1.pdf




####model coding
##random effect structure is NESTED:
###(1|Colony/Microcolony/Time)
#This means a random intercept for each timepoint....
#nested within microcolony...
#nested within colony

#Two possible packages to use:
#1. nlme package, function lme
#2. lme4 package, function lmer 
#I will code use lme in nlme;
#only nlme::lme gives p-values for the fixed effects coefficients
#lme4 package authors philosophically opposed to p values

library (nlme)
M0<-lme(PRate ~ Infection * Pollen * Time , data=Data, #data frame to use
        random= ~1|Colony/Microcolony/Time) 
summary(M0)
#here is the same model with lme4 
#library(lme4)
#AltModel<-lmer(PRate ~ Infection * Pollen * Time + (1|Colony/Microcolony/Time), data=Data)

#exclude 3-way interaction
M1<-update (M0, ~. - Infection:Pollen:Time)
anova(M0, M1) #oops- doesn't work because model fit with "REML" not "ML"
library(bbmle)
AICtab(M0, M1) #M1 better, has lower AIC
summary(M1) 
#drop Infection:Time
M2<- update(M1,  ~. - Infection:Time)
AICtab(M0, M1, M2) #M2 is best, keep simplifying
summary(M2)
#drop Pollen:Time 
M3<- update(M2,  ~. - Pollen:Time)
AICtab(M3, M2) #M3 better, keep simplifying
summary(M3)
# drop Infection:Pollen
M4<-  update(M3,  ~. -Infection:Pollen  )
AICtab(M4, M3) #M4 better, keep simplifying
summary(M4)
#None significant-- infection, pollen , time
M5<-update(M4, ~. -Time )

library(car)
Anova(M5)

PollenModel<-M5
M5
summary(M5) #this is more complete for the fixed effects coefficients and std errors
#M5<-lme(PRate ~ Infection + Pollen , data=Data, random= ~1|Colony/Microcolony/Time) 
####plot model results
#don't think we need this so much if no effect of the main treatments
lsmeans(M5, ~Pollen)
#pretty big standard errors, confidence intervals spanning 0 in both cases
#we probably need a better way to measure pollen consumption 
#to increase signal to noise ratio



#http://glmm.wikidot.com/faq
#use this script to generate a reference grid and plot confidence intervals

