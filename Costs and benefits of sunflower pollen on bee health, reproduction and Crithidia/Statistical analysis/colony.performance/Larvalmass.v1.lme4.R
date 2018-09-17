####2016.03#########
#analysis of sunflower microcolony performance

#####Larval production################
####Source data: use summary created by LS Adler
#EPY added back colonies from which new adults emerged -- 
#colonies 3, 11, 23, 25, 33
#coded as "1" for pupal production

#rm(list=ls()) #clear memory

setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/jess stats/colony.performance")

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

###HERE I WILL DEAL WITH THE MICROCOLONIES THAT MADE LARVAE
Makers
length(Makers$Microcolony) #41 microcolonies
Boxplot(meanlarvwt~Infection*Pollen, Makers)
#no way that this meets assumptions of equal variance
TABLE<-ddply(Makers,~Infection*Pollen, 
               summarise,numlarvae=mean(larvnum),
               meanmass=mean(larvwt, na.rm=TRUE),
             samplesize=count(larvnum))
xtabs( ~Pollen, Makers)
xtabs( ~Infection, Makers)


#plenty of zeroes
#check distribution
qqp(Perf03$meanlarvwt, "norm") #OK but does not account for different variance!
#anyway the treatment effect is dumb clear in the boxplot,
#so I will just run the model anyway because we have to analyze somehow


#try lognormal
# lnorm means lognormal
qqp(Perf03$meanlarvwt, "lnorm") #we don't need such broad bands
#log-normal would be a biologically appropriate distribution

Buck<-subset(Perf03, Pollen=="B")
Sun<-subset(Perf03, Pollen=="S")

#try gamma
gamma <- fitdistr(Makers$meanlarvwt, "gamma")
gamma <- fitdistr(Makers$meanlarvwt, "gamma")
qqp(Makers$meanlarvwt, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])
#like this b/c allows variance to increase w/ mean


#sample gamma fit
#(m <- glmer(y ~ x + (1 | id), family = Gamma))

#first try poisson model:
M1<- glmmTMB(meanlarvwt ~ Infection*Pollen  + (1|Colony), 
                    family = "Gamma",
                    data=Makers)
#Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')
# This often occurs when a model is overparameterized

#Try in lme4 instead
M1<- glmer(meanlarvwt ~ Infection*Pollen  + (1|Colony), 
             family = "Gamma",
             data=Makers)

summary(M1)
drop1(M1, test = "Chisq")
# Df     AIC     LRT Pr(Chi)
# <none>              -196.13                
# Infection:Pollen  1 -197.99 0.14168  0.7066
Anova(M1)

#drop interaction term
M2<-update(M1, ~. - Infection:Pollen)
summary(M2)
Anova(M2)
# Chisq Df Pr(>Chisq)    
# Infection  0.0176  1  0.8944439    
# Pollen    11.4439  1  0.0007173 ***

#good-- big effects of pollen


################### Plots #############
Lsm<-lsmeans(M2, ~Pollen)
Pairs<-contrast(Lsm, "pairwise") 
Pairs
# contrast estimate       SE df z.ratio p.value
# B - S    66.76646 19.73657 NA   3.383  0.0007

##Sunflowers distinguish themselves

##Get letters to plot
Lets<-cld(Lsm, Letters = letters) 
Lets
Lets$.group<-gsub(" ", "", Lets$.group) ##remove spaces

################## PLOT ON SCALE OF RESPONSE ##########
#Let's invert to original mean+SE and mean-SE
Lsmeans.df<-as.data.frame(summary(lsmeans(M2, ~Pollen)))
Lsmeans.df$meanlarvwt<-1/(Lsmeans.df$lsmean)
Lsmeans.df$btlo<-1/(Lsmeans.df$lsmean + Lsmeans.df$SE) #lower bound-- note inverse link
Lsmeans.df$bthi<-1/(Lsmeans.df$lsmean - Lsmeans.df$SE) #upper bound-- note inverse link

Lsmeans.df

#Merge with letters
Lets2<- dplyr::select(Lets, Pollen, .group)
Lets2

Lsmeans.df.merged<- join(Lsmeans.df, Lets2, by = "Pollen")

Lsmeans.df.merged
# Pollen   lsmean        SE df asymp.LCL asymp.UCL meanlarvwt       btlo        bthi .group
# 1      B 87.41157 19.340245 NA  49.50539 125.31775 0.01144013 0.01469047 0.009367522      b
# 2      S 20.64511  4.509244 NA  11.80716  29.48307 0.04843761 0.06197373 0.039754542      a
#More than double larval mass

#Manually switch the letters to avoid inversion
Lsmeans.df.merged$.group[Lsmeans.df.merged$Pollen=="B"]<-"a"
Lsmeans.df.merged$.group[Lsmeans.df.merged$Pollen=="S"]<-"b"


Lsmeans.df<-Lsmeans.df.merged #overwrites old

#x-axis:
pollen.axis<-c("Buck", "Sun")
#y-axis:
ylabel<-"Mean larval mass (g)"

#ready to plot?
p<-ggplot(Lsmeans.df, aes(x = Pollen, y = meanlarvwt)) + 
  #Option to sprinkle raw data... 
  geom_point(data = Makers, aes( x = Pollen, y = meanlarvwt), position = "jitter", 
             color = "blue", shape = 25, size = 0.5, alpha = 0.4)+
  geom_pointrange(aes(ymin = btlo, ymax = bthi), size = 1, color = "black") +
  scale_x_discrete(labels=pollen.axis) +
  ylab(ylabel) + #y label
  xlab("Pollen")
p
#Add text
P.lmass<- p + geom_text(aes(y = bthi + 0.4 * max(bthi), label = .group),
                     size = 8)
P.lmass

##Export
# setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/jess stats/colony.performance")
# ggsave("LARVALMASS.v1.lme4.pdf", height = 5, width =5)




