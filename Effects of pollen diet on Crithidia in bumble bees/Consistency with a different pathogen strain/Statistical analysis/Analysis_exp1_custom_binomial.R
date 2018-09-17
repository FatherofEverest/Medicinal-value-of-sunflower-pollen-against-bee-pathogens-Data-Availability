####analysis of POLLEN EXPERIMENT_1_JC Roulston Crithidia Strain ####### 
###individual bees.... JONATHAN GIACOMINI.... 2 POLLEN (Sun + Buck ) DIETS EFFECTS ON CRITHIDIA COUNTS#
#JG:
data1<-read.csv("/Users/rg9403/Dropbox/Sunflower Pollen Experiments/Exp1/Sunflower Pollen Exp 1_Analysis 1 data.csv",header=TRUE)
#epy/Google Drive
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/NCSU_Crithidia/Exp1_JCRoulstonStrain_BuckSun")
data1<-read.csv("Sunflower Pollen Exp 1_Analysis 1 data.csv")
str(data1) 
#following should be factors:
data1$BeeID<-as.factor(data1$BeeID)
data1$Date<-as.factor(data1$Date)

#scout data
hist(data1$WingSize) #good, no pterodactyls
hist(data1$Count)

Zeroes<-subset(data1, Count==0)
length(Zeroes$BeeID)
#ouch, 108 zeroes!
length(data1$BeeID)
108/148
#>70% zeroes!
#if you toss the zeroes, you only have 40 bees left across 
#2 treatments *3 colonies * 2 inoculation dates=12 combinations!

#This is going to be touch with anything except a binomial model
#However, seeing as when sunflower counts>0 they are very low,
#it might be most meaningful to define infection YES at some threshold,
#say count>3 Crithidia cells in 0.02uL gut extract

xtabs(~Treatment, data=Zeroes)
#twice as many sunflower had zeroes
#Treatment
#BUCK  SUN 
#38   70 
xtabs(~ColonyID, data=Zeroes)
#ColonyID
#JGA JGB JGC 
#39  40  29 
#all pretty even, zeroes not driven by one colony
xtabs(~Date, data=Zeroes) #again pretty even,
#Date
#15336 15338 
#57    51 
#so not as if you just had a wacky inoculation mix on one of the days

#let's plot the distribution
library(MASS)
library(car)
#scout distribution, see helpful script by julia pilowsky
#http://www.juliapilowsky.com/mixedmodels/
Data1<-data1
qqp(Data1$Count, "norm") ##forget it
# lnorm means lognormal
qqp(Data1$Count, "lnorm") #forget it
#Poisson:
poisson <- fitdistr(Data1$Count, "Poisson")
qqp(Data1$Count, "pois", poisson$estimate) #forget it
#not good, way over-dispersed

#negative binomial
nbinom <- fitdistr(Data1$Count, "Negative Binomial") #can't even fit 
qqp(Data1$Count, "nbinom", size = nbinom$estimate[[1]],
    mu = nbinom$estimate[[2]])
#looks like it only used the 40 bees with count>0???
#try with N+1 transformation (really a translation)
Data1$Count.t<-Data1$Count +1
nbinom <- fitdistr(Data1$Count.t, "Negative Binomial") #can't even fit 
qqp(Data1$Count, "nbinom", size = nbinom$estimate[[1]],
    mu = nbinom$estimate[[2]]) #ooftah

Boxplot(Data1$Count~Data1$Treatment) #median 0 in both groups

#define new columns called "Infection" with count>3
#mydata$NewTemp <- ifelse(mydata$Temp>0, 1, 0)
data1$Infection <- ifelse(data1$Count>3, 1, 0) #sets infection threshold
#alternatively:
#mydf$TempBin <- as.numeric(mydf$Temp > 70)
data1$Infectionbin <- as.numeric(data1$Count > 3)
Compare<-cbind(data1$Infection, data1$Infectionbin)
View(Compare)

#Now run binary model
Data1<-na.omit(data1) #The following models don't work unless I omit the NAs
library(glmmADMB)
library(lme4)
Data1$Infection
str(Data1$ColonyID)
str(Data1$WingSize)
Fit1<-glmmadmb(Infection ~ Treatment  + WingSize +(1|ColonyID)+(1|Date), 
               family = "binomial",
               data=Data1)
summary(Fit1)
#Log-likelihood: -54.5819 
plot(Fit1$residuals)
#drop wing size
Fit2<-update(Fit1, ~. -WingSize)
summary(Fit2)
library(bbmle)
AICtab(Fit1, Fit2) #prefer Fit2

Anova(Fit2) #effect of treatment

library(lsmeans)
My.lsmeans<-lsmeans(Fit2, ~Treatment)
LSdf<-summary(My.lsmeans)
#manual transformation using inverse logit
My.lsmeans
library(boot)
LSdf$btmean<-inv.logit (LSdf$lsmean)
LSdf$btmean #these are probabilities
LSdf$btSElo<-inv.logit(LSdf$lsmean - LSdf$SE)
LSdf$btSEhi<-inv.logit(LSdf$lsmean + LSdf$SE)
library(ggplot2)
P0<-ggplot(LSdf, aes(x=Treatment, y=btmean)) + 
  geom_bar(position=position_dodge(), stat="identity", size=2,colour="black" ) +
  geom_errorbar(aes(ymin=btSElo, ymax=btSEhi), size=2, width=0.2)+
  theme_classic(base_size=25)
P0
pollen.axis<-c( "Buckwheat", "Sunflower")
ylabel<-expression(bold(Probability~of~infection))
FigTF<-P0 + theme(text=element_text(face="bold") , line=element_line(size=2)) +
   ylab(ylabel) + xlab("Pollen type") + 
  scale_x_discrete( labels=pollen.axis) +
  scale_y_continuous(expand = c(0,0))
pdf("custom_binomial_v0_manual_standard.error.pdf")
FigTF
dev.off() ###THIS HAS MANUALLY COMPUTED SE on probability scale 
###I would prefer to show this rather than the below

#transform to scale of probabilities using lsmeans
My.Lsmeans2<-lsmeans(Fit2, ~Treatment, type="response")
My.Lsmeans2
#My.Lsmeans2
#Treatment       prob         SE df   asymp.LCL asymp.UCL
#BUCK      0.33069021 0.09239153 NA 0.178988911 0.5282406
#SUN       0.02324546 0.02440700 NA 0.002885978 0.1636592
#only 2% of sunflower bees expected to have count>3
library(ggplot2)
Lsmeans.df<-summary(My.Lsmeans2)
str(Lsmeans.df)
View(Lsmeans.df)


###############Plotting################
pollen.axis<-c( "Buckwheat", "Sunflower")

plot0<- ggplot(Lsmeans.df, aes(x=Treatment, y=prob)) + 
  geom_bar(position=position_dodge(), stat="identity", colour="black" ) +
  geom_errorbar(aes(ymin=prob-SE, ymax=prob+SE),
                #####note you can also
                #geom_errorbar(aes(ymin=asymp.LCL, ymax=asymp.UCL),
                size=2, width=0.2,  # Width of the error bars
                position=position_dodge(0.9))
plot0 # Looks ok, let's refine it now
plot01<-plot0 + theme_bw() #remove colored background
#ylabel<- expression(bold(paste("Parasite load (ln (cells * 0.02 ",mu, L^-1, "))", sep="") )) #check for encoding
ylabel<-expression(bold(Probability~of~infection))
plot02<-plot01 + ylab(ylabel) #y label
plot02
#relabel x axis:
plot05<-plot02 + xlab("Pollen type") + # x axis label
  scale_x_discrete(breaks=levels(Lsmeans.df$Treatment) , labels=pollen.axis)
plot05
#Larger font:
plot06<-plot05 + theme(text=element_text (size=25, face="bold") )
plot06
#Remove gridlines
neater<- plot06 +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.line = element_line(size = 1.5, colour = "black"))
neater
#axis scaling
Zeroed<- neater + coord_cartesian(ylim = c(0, 0.43))
Zeroed
#if we want y axis to start at 0:
#this removes space between x axis and bottom of bar (use for non-log scale)
Justified<- neater +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.43))
Justified
fig1<-Zeroed
fig1 #here the bars show standard error of probability
pdf("custom_binomial_v0_standard.error.pdf", height=7, width=7, paper="special")
fig1
dev.off()







#Here's JG's original negative binomial model
#but see the plots above-- negative binomial
#distribution fits response var. horridly!
############## Negative Binomial Fit ##################
data1$Date<-factor(data1$Date) # Change Date to a factor
Data1<-na.omit(data1) #The following models don't work unless I omit the NAs

Fit1<-glmmadmb(Count ~ Treatment  + WingSize +(1|ColonyID)+(1|Date), 
               family = "nbinom",
               data=Data1)
summary(Fit1)
library(car)
Anova(Fit1)
# Analysis of Deviance Table (Type II tests)
# 
# Response: Count
# Df   Chisq Pr(>Chisq)    
# Treatment   1 30.1764  3.945e-08 ***
#   WingSize    1  9.1121   0.002539 ** 
#   Residuals 141

# We can knockout Date and colony ID and then compare models
Fit2<-update(Mod1, ~. - Date)
summary(Fit2)
Anova(Fit2)

Fit3<-update(Mod1, ~. - ColonyID)
summary(Fit3)
Anova(Fit3)


AIC(Fit1) #459.914 #### Best fit...Let's keep Both Colony ID and Date as random effects
AIC(Fit2) #470.988
AIC(Fit3) #470.988

######## PLOT ######## Not so good - the Lsmean for sun is negative 
library(lsmeans)
My.Lsmeans2 <- lsmeans(Fit1, ~Treatment) 
My.Lsmeans2
library(ggplot2)
Lsmeans.df<-summary(My.Lsmeans2)
str(Lsmeans.df)
View(Lsmeans.df)
#label for x axis
pollen.axis<-c( "Buckwheat", "Sunflower")

plot0<- ggplot(Lsmeans.df, aes(x=Treatment, y=lsmean)) + 
  geom_bar(position=position_dodge(), stat="identity", colour="black" ) +
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE),
                size=1, width=0.2,  # Width of the error bars
                position=position_dodge(0.9))
plot0 # Looks ok, let's refine it now
plot01<-plot0 + theme_bw() #remove colored background
ylabel<- expression(bold(paste("Parasite load (ln (cells * 0.02 ",mu, L^-1, "))", sep="") )) #check for encoding
plot02<-plot01 + ylab(ylabel) #y label
plot02
#relabel x axis:
plot05<-plot02 + xlab("Pollen type") + # x axis label
  scale_x_discrete(breaks=levels(Lsmeans.df$Treatment) , labels=pollen.axis)
plot05
#Larger font:
plot06<-plot05 + theme(text=element_text (size=25, face="bold") )
plot06
#Remove gridlines
neater<- plot06 +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.line = element_line(size = 1.5, colour = "black"))
neater
#axis scaling
Zeroed<- neater + coord_cartesian(ylim = c(0, 6.5))
Zeroed
#if we want y axis to start at 0:
#this removes space between x axis and bottom of bar (use for non-log scale)
Justified<- neater +
  scale_y_continuous(expand = c(0,0), limits = c(0,6.5))
Justified
fig1<-Justified
fig1



