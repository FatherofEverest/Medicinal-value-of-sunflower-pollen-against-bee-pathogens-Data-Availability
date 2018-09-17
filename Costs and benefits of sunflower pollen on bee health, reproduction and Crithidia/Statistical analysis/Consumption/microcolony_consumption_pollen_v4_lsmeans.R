#update 170605: larger symbols
#update 160312: use lsmeans instead of effects for plot
#update 2016.03.01: time as fixed effect only

#load packages
library(nlme)
library(lme4)
library(car)
library(plyr)
library(lsmeans)
library(ggplot2)

rm(list=ls())
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/jess stats/Consumption")
Data<-read.csv(
  "C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony//Pollen consumption trials/updated corrected consumption values only.csv",
  na.strings = ".")
#View(Data)
str(Data)
Data$Colony<-Data$Start.Day
Data$Time<-Data$Day
Data$Infection<-Data$Infect.
#change reference level
library(plyr)
Data$Infection <- relevel(Data$Infection, ref = "U")

Data$Microcolony<-as.factor(Data$Micro...)
str(Data$Microcolony)
#eliminate colonies 12 and 16, these were uninfected treatments that had infection.
#also eliminate colonies 50 and 54-- died within first week
Data0<-subset(Data, !Microcolony=="12" & !Microcolony=="16" & !Microcolony=="50" & !Microcolony=="54")
str(Data0)
Data0$Microcolony
Data<-droplevels(Data0)
Data$MicrocolonySize<-Data$X..of.bees
Data$NRate<-Data$Nectar.per.bee.per.hour
Data$PRate<-Data$Pollen.per.bee.per.hour
str(Data)

###add size dimorphism (optional)
#I prefer not to do have this in model because there are so many values missing
SizeData<-read.csv("C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/Pollen consumption trials/updated survival.csv")#rename some columns
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
str(Data$Microcolony)
str(Data$Colony)

#####data inspection (plotting)
#see some superb examples of plotting by Ben Bolker
#look for his "Lab 1 mixed models"
#http://ms.mcmaster.ca/~bolker/classes/uqam/mixedlab1.pdf
hist(Data$PRate)
library(car)
Boxplot(Data$PRate)
#get rid of observation 864, 688-- appears to be data entry error
Data<-Data [-c(826, 657),]
Boxplot(Data$PRate) #more reasonable

hist((Data$PRate)) #GOOD -- distribution much improved
hist(log(Data$PRate)) #won't work now that I have brought back the negative

Boxplot(Data$PRate~Data$Pollen) 
#here it definitely looks like sunflower was more consumed

Data<-droplevels(Data)
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
M0randomT<-lme(PRate ~ Infection * Pollen * Time , data=Data, #data frame to use
        random= ~1|Colony/Microcolony/Time, method="ML") 
summary(M0randomT)
Anova(M0randomT)

#Set contrasts to test interaction term
contrasts <-c("contr.sum", "contr.poly")
M0<-lme(PRate ~ Infection * Pollen * Time , data=Data, #data frame to use
               random= ~1|Colony/Microcolony, method="ML") 
summary(M0)
Anova(M0)
library(bbmle)
AICtab(M0, M0randomT) #differ by 2 aic units; use simpler model
anova(M0, M0randomT) #no difference

#start simplification from here
summary(M0)
#exclude Infection:Time
M01<-update(M0, ~. - Infection:Time)
anova(M01, M0)


logLik((M0))
logLik(M01) #why are these exactly the same?
AIC(M0)
AIC(M01)
anova(M0, M01)
summary(M01)
Anova(M01)



#3 way interaction:
M02<-update(M01, ~. - Infection:Pollen:Time)
anova(M02, M01) #p=0.10 in likelihood ratio test

Anova(M02)

#likelihood ratio test for Infection:Pollen
M03<-update(M02, ~. - Infection:Pollen)
anova(M02, M03) #p=0.06, drop
Anova(M03)
Anova(M03, type = 3)
summary(M03)
#M03<-lme (PRate ~ Infection + Pollen + Pollen:Time , data=Data, random= ~1|Colony/Microcolony, method="ML") 



#here is the same model with lme4 
Data<-droplevels(Data)
lmer03<- lmer(PRate ~ Infection + Pollen + Pollen:Time + (1|Colony/Microcolony), data=Data)




#and starting from scratch in lme4
library(lme4)
AltModel<-lmer(PRate ~ Infection * Pollen * Time + (1|Colony/Microcolony), data=Data)
Anova(AltModel)

Alt2<- lmer(PRate ~ Infection * Pollen * Time + Colony + (Time|Microcolony), data=Data, REML=FALSE)

Alt2a<- lmer(PRate ~ Infection * Pollen * Time + (1|Colony) + (Time|Microcolony), data=Data, REML=FALSE)
#above does not account for nesting of microcolony
Anova(Alt2a) #see the same fixed effects as in the lme model (M03)
AICtab(Alt2a, M0) #lmer model with random slopes has (much) lower AIC
#however, we have time as a fixed effect already...
anova(Alt2, Alt2a) #higher log-likelihood with colony as random (Alt2a)
AICtab(Alt2, Alt2a) #lower AIC with colony as fixed effect (Alt2)
#Lynn wanted colony as random (philosophical)
#try having microcolony within colony and random slopes for mc by time
Alt2b<- lmer(PRate ~ Infection * Pollen * Time + (1|Colony/Microcolony) + (Time|Microcolony), data=Data, REML=FALSE)
#warning message, model failed to converge
summary(Alt2b) #"unable to evaluate...."
Anova(Alt2b) #car package is very cooperative:)

Alt2b<- lmer(PRate ~ Infection * Pollen * Time + (1|Colony/Microcolony) + (Time|Microcolony), data=Data, REML=TRUE)
#"Model is nearly unidentifiable: large eigenvalue ratio"-- probably too many random effects
Alt2b<- lmer(PRate ~ Infection * Pollen * Time + (1|Colony/Microcolony) + (Time|Microcolony), data=Data, REML=TRUE)


#try flipping order of mc and time
Alt2c<- lmer(PRate ~ Infection * Pollen * Time + (1|Colony/Microcolony) + (Microcolony|Time), data=Data, REML=FALSE)
#outright error, too many random effects


######Here is simplification using colony as fixed effect#####
summary(Alt2)
Anova(Alt2)

#exclude 3-way interaction
M1<-update (Alt2, ~. - Infection:Pollen:Time)
anova(Alt2, M1) #only works if REML=FALSE for both models
library(bbmle)
AICtab(Alt2, M1) #M1 better, has lower AIC, difference miniscule
summary(M1) 
Anova(M1)
#drop Infection:Time
M2<- update(M1,  ~. - Infection:Time)
AICtab(M0, M1, M2) #M2 is best, keep simplifying
summary(M2)
Anova(M2)
#drop Infection:Pollen
M3<- update(M2,  ~. - Infection:Pollen)
AICtab(M3, M2) #M3 better, keep simplifying
summary(M3)
Anova(M3)
#Good, we end up with about the same result with lme,
#except now we have colony among the fixed effects
#highly significant effect of colony, by the way

summary(M3) #this is more complete for the fixed effects coefficients and std errors
#M5<-lme(PRate ~ Infection + Pollen , data=Data, random= ~1|Colony/Microcolony/Time) 
####plot model results
#don't think we need this so much if no effect of the main treatments
lsmeans(M3, ~Pollen*Time)

AICtab(M3, M03) #the lmer model M3 scores much better here
logLik(M3)
logLik(M03) #again higher scores for the lmer "M3" model with random slopes

####Now do the same lmer simplification using colony as random effect######
Alt2a<- lmer(PRate ~ Infection * Pollen * Time + (1|Colony) + (Time|Microcolony), data=Data, REML=FALSE)
Anova(Alt2a)
#exclude 3-way interaction
M1a<-update (Alt2a, ~. - Infection:Pollen:Time)
anova(Alt2a, M1a) #only works if REML=FALSE for both models
library(bbmle)
AICtab(Alt2a, M1a) #M1 better, has lower AIC, difference miniscule
summary(M1a) 
Anova(M1a)
#drop Infection:Time
M2a<- update(M1a,  ~. - Infection:Time)
AICtab(M1a, M2a) #M2 is best, keep simplifying
summary(M2a)
Anova(M2a)
#drop Infection:Pollen
M3a<- update(M2a,  ~. - Infection:Pollen)
AICtab(M3a, M2a) #M3 better, keep simplifying
summary(M3a)
Anova(M3a)

logLik((M3))
logLik(M3a) #slightly worse for M3a with colony as random
AICtab(M03, M3, M3a) #again worse for M3a with colony as random
#random slope model is much better in terms of AIC, but most df

#reasons to go with the "lme" fit  (M03) because 
#despite higher AIC, 
#that uses the same random effect structure as the counts and survival
#also we already have time as a fixed effect, so don't want to have time as random effect also
MyModel<- M03


#first with "lsmeans"
timevec=seq(0,42,1)
timevec
library(lsmeans)
ls2plot<- summary(lsmeans(MyModel, ~Pollen*Time, at=list(Time=timevec)))
#View(ls2plot)
#change scale to ug/bee/h
ls2plot$meanug<-ls2plot$lsmean*10^6
ls2plot$SEug<-ls2plot$SE*10^6

library(cowplot)
levels(ls2plot$Pollen)
p0<-ggplot(data=ls2plot, aes(y=meanug, x=Time, lty=Pollen))  +
  geom_line(stat="identity", size=2)+ 
  #geom_errorbar(limits, position=dodge) +
  geom_ribbon(aes(ymin=meanug - SEug, ymax=meanug + SEug), alpha=0.5) #pretty!
ylab<-expression(bold(Pollen~consumption~(mu*g~bee^-1*h^-1)))
p1<-p0 +  theme_cowplot() +
  xlab("Days since inoculation")+ylab(ylab)+
  scale_linetype_discrete(name="Pollen", labels=c("Buckwheat", "Sunflower")) +
    theme(text=element_text(face="bold"), 
        line=element_line(size=2)) +
  theme(axis.line.x=element_line(size=2),
  axis.line.y=element_line(size=2))
  p1
 p2<- p1 + theme(legend.key.size=unit(1.5, "cm"))
p3<- p2 + scale_x_continuous(expand = c(0,0))#removes gap



setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/jess stats/Consumption")
# pdf("Pollen_consumption_v4_lsmean.pdf", width=7, height=7, paper="special")
# p3
# dev.off()

##note that we lose the ribbons with this
#postscript("Pollen_consumption_v4_lsmean.eps", horizontal = FALSE, onefile = FALSE, 
 #          paper = "special", height=7,width=7)
p3
#ggsave("Pollen_consumption_v4_lsmean.eps", device=cairo_ps, height=7,width=7)

#dev.off()
 
#option to overlay raw mean +/- SE
levels(Data$Pollen)
levels(ls2plot$Pollen)
p3 + geom_point(data=Data, aes(x=Time, y=10^6*PRate))
#a little too swamped
library(plyr)
Data$Timev
cdata <- ddply(Data, c("Pollen", "Time"), summarise,
               N    = sum(!is.na(PRate)),
               mean = 10^6*mean(PRate, na.rm=TRUE),
               sd   = 10^6*sd(PRate, na.rm=TRUE),
               se   = 10^6*sd / sqrt(N)
)
#View(cdata)
p4<-p3 + geom_point(data=cdata, 
    aes(x=Time, y=mean, shape = Pollen))+
  xlim(c(0,41)) + ylim(c(-200,900))+
  scale_shape_discrete(name=NULL, labels=c("Buckwheat", "Sunflower"))
p4
#ggsave("Pollen_consumption_v5_cowplot.pdf", height = 5, width = 8)

#Export to main figure directory
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Manuscript Drafts/Figures")
#ggsave("Pollen_consumption_v5_cowplot.pdf", height = 5, width = 8)


