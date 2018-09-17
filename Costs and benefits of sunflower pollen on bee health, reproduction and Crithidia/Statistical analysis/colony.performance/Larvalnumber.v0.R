####2016.03#########
#analysis of sunflower microcolony performance

#####Larval production################
####Source data: use summary created by LS Adler
#EPY added back colonies from which new adults emerged -- 
#colonies 3, 11, 23, 25, 33
#coded as "1" for pupal production

rm(list=ls()) #clear memory

library(dplyr)
library(car)
library(nlme)
library(lme4)
library(glmmADMB) #takes some work to install this
library(bbmle)
library(ggplot2)

Perf0<- read.csv('C:/Users/Evan/Dropbox/Jess Sum15/jess stats/colony.performance/performance_sum_SAS_CALLOWS.ADDED.csv',
                 na.strings=".")
Bees<-read.csv("C:/Users/Evan/Dropbox/Jess Sum15/Pollen consumption trials/updated survival.csv")
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
#too many zeroes, poor fit

#try gamma
gamma <- fitdistr(Perf03$larvnum, "gamma")
gamma <- fitdistr(1+Perf03$larvnum, "gamma")
qqp(1+ Perf03$larvnum, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])
#OK, but this gamma is for continuous distribution

#let's compare Poisson and negative binomial models
library(glmmADMB)
#first try poisson model:
fitPoiss<- glmmadmb(larvnum ~ Infection*Pollen  + (1|Colony), 
                    family = "poisson",
                    data=Perf03)
fitNB<- glmmadmb(larvnum ~ Infection*Pollen  + (1|Colony), 
                    family = "nbinom",
                    data=Perf03)
library(bbmle)
plot(fitPoiss$residuals)
plot(fitNB$residuals) #not much difference
AICtab(fitPoiss, fitNB)
#HUGE edge to nbinom
logLik(fitPoiss)
logLik(fitNB) #NB MUCH BETTER

summary(fitNB)
Anova(fitNB)
#drop interaction
NB1<-update(fitNB, ~. - Infection:Pollen)
Anova(NB1)

#Let's plot this
library(lsmeans)
My.lsmeans<-as.data.frame(summary(lsmeans(NB1, ~Pollen)))
My.lsmeans
margmeans2<-as.data.frame(summary(lsmeans(NB1, ~Pollen)), type="response")
margmeans2
?lsmeans #type="response" is supposed to exponentiate, but no luck
LCOUNTS<-ddply(Perf03,~Infection*Pollen, 
               summarise,numlarvae=mean(larvnum),
               meanmass=mean(larvwt, na.rm=TRUE))
LCOUNTS #these are the arithmetic means, clearly much higher 
#than in the lsmeans output
#I will exponentiate manually
margmeans<-My.lsmeans
margmeans$btmean<-exp(margmeans$lsmean)
margmeans$rawhi<-margmeans$lsmean + margmeans$SE
margmeans$rawlo<-margmeans$lsmean - margmeans$SE
margmeans$bthi<-exp(margmeans$rawhi)
margmeans$btlo<-exp(margmeans$rawlo)
##now the plot
library(ggplot2)
#make a data frame from the lsmeans object
Lsmeans.df<-margmeans
str(Lsmeans.df)
View(Lsmeans.df)
#label for x axis
pollen.axis<-c( "Buckwheat", "Sunflower")

plot0<- ggplot(Lsmeans.df, aes(x=Pollen, y=btmean)) + 
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", size=2 ) +
  geom_errorbar(aes(ymin=btlo, ymax=bthi),
                size=2, width=0.2,  # Width of the error bars
                position=position_dodge(0.9))+
  theme_classic(base_size=24)
#where are my axes?
p1<-plot0 + ylab("Number of Larvae") + xlab("Pollen")+
   theme (axis.line.y=element_line(color="black", size=2),
           axis.line.x=element_line(color="black", size=2)) +
  theme(text=element_text(face="bold"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 7.9)) 
p1
p2<-p1+ scale_x_discrete(labels=c("Buckwheat", "Sunflower"))
p2
fig1<-p2
#export to PDF
#dev.new() #sometimes this helps with errors
setwd("C:/Users/Evan/Dropbox/Jess Sum15/jess stats/colony.performance")
pdf("LARVALNUMBER.pdf",height=6, width=5, paper='special') 
fig1
dev.off()

#for eps
setEPS()
postscript("LARVALNUMBER.eps", horizontal = FALSE, onefile = FALSE, 
           paper = "special", height = 6, width = 5)
fig1
dev.off()

