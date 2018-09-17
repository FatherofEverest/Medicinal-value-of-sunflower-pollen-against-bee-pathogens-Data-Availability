####2016.03#########
#analysis of sunflower microcolony performance

#####Pupal production################
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
View(Perf02)
#rename(Perf02, Colony=sourcecol) #not working
Perf02$Colony<-Perf02$sourcecol
str(Perf02) #duplicated, renamed the column
Perf03<-Perf02
str(Perf03$infect) #good , 2 levels
Perf03$Infection<-Perf03$infect
Perf03$Infection<-relevel(Perf03$Infection, ref="U")
str(Perf03) #renamed the column

Perf03$Pollen<-Perf03$pollen
?rename
#check pupae
Perf03$pup.binary
length(Perf03$pup.binary)
library(car)
library(plyr)
?summarize #to summarize by a grouping variable
avg.probs<-ddply(Perf03,~Pollen, summarise,Production=mean(pup.binary))
avg.probs #28% of sunflower colonies made pupae; NO buckwheat colonies made pupae!
xtabs(pup.binary~Pollen, data=Perf03)
#11 sunflower microcolonies made pupae, 0 buckwheat m-colonies
xtabs(pup.binary~Infection, data=Perf03)
#5 infected, 6 uninfected made pupae
xtabs(~Pollen, data=Perf03)

xtabs(pup.binary~Colony, data=Perf03)
#0 in JML-3, 1 in JML-4, 5 each in JML-1 and JML-2

#ready to run binomial model with Colony as random effect?
GHQ <- glmer(pup.binary ~ Pollen*Infection +  (1 | Colony), data = Perf03,
             family = binomial(link = "logit"), nAGQ = 100)  # Set nAGQ to # of desired iterations
#warning maybe b/c buckwheat has no pupae
summary(GHQ)
Anova(GHQ)

G1<-update(GHQ, ~. - Pollen:Infection)
#convergence error
summary(G1)
Anova(G1) #not happy, "variance-covariance matrix computed 
#from finite-difference Hessian is
#not positive definite or contains NA values"

G1default<-glmer(pup.binary ~ Pollen+Infection +  (1 | Colony), data = Perf03,
                 family = binomial)  
summary(G1default)

G2<-glmer(pup.binary ~ Pollen+Infection +  (1 | Colony), data = Perf03,
                 family = binomial(link = "logit"), nAGQ = 100)  # Set nAGQ to # of desired iterations

#try in glmmADMB
AD1<-glmmadmb(pup.binary ~ Pollen+Infection +  (1 | Colony), data = Perf03,
          family = "binomial")  # Set nAGQ to # of desired iterations
summary(AD1)
Anova(AD1)
AD2<-glmmadmb(pup.binary ~ Pollen+Infection +  (1  ), data = Perf03,
              family = "binomial")  # Set nAGQ to # of desired iterations
summary(AD2)
Anova(AD2)


#What about removing colony? This would be a fixed effects model
Norand<-glm(pup.binary ~ Pollen+Infection + Colony , data = Perf03,
             family = binomial(link = "logit"))
#Warning- glm.fit: fitted probabilities numerically 0 or 1 occurred 
#model is having trouble with zero successes among buckwheat
summary(Norand) #non-significant coefficients... BUT
Anova(Norand) #????? now the model is confident 
#that colony and pollen are important 
#epy: I don't think the model is working--
#note that Beta for JML-3 is -1.96 (NO PUPAE), 
#but -2.2 for JML-4 (1 PUPA) 
#Shouldn't JML-3 have the most negative estimate??

#give up, what about just a chi-square test?
chisq.test(x=Perf03$Pollen, y=Perf03$pup.binary, simulate.p.value=TRUE)
chisq.test(x=Perf03$Infection, y=Perf03$pup.binary, simulate.p.value=TRUE)
chisq.test(x=Perf03$Colony, y=Perf03$pup.binary, simulate.p.value = TRUE)





#We could also use just the subset of Colonies JML1:2
#THIS APPROACH ENDORSED BY LSA ON 2016.03.02####

Producers<-subset(Perf03, !Colony=="JML-3" & !Colony=="JML-4")
Producers$Microcolony

#Run the chi-square tests again on the subset
chisq.test(x=Producers$Pollen, y=Producers$pup.binary, simulate.p.value=TRUE)
chisq.test(x=Producers$Infection, y=Producers$pup.binary, simulate.p.value=TRUE)
chisq.test(x=Producers$Colony, y=Producers$pup.binary, simulate.p.value = TRUE)

#generalized linear model (no random effects) on the JML 1 & 2 subset
Producers<-droplevels(Producers)
GLM.12<-glm(pup.binary ~ Pollen+Infection  , data = Producers,
            family = binomial(link = "logit"))
summary(GLM.12)
Anova(GLM.12)


###plot results
library(lsmeans)
lsmeans(GLM.12, ~Pollen, type="response")
#this looks a little wacky, maybe just plot the values...
#let's get an idea with some tables:
xtabs(pup.binary~Pollen, data=Producers)
xtabs(pup.binary~Infection, data=Producers)
xtabs(pup.binary~Colony, data=Producers)
xtabs(pup.binary~Pollen, data=Perf03)
xtabs(pup.binary~Infection, data=Perf03)
xtabs(pup.binary~Colony, data=Perf03)
xtabs(pup.binary~Colony, data=Producers)
xtabs(~Pollen, data=Producers)

#now plot just the raw proportions
#repeat chi-square tests with micro'
library(ggplot2)
avg.probs<-ddply(Producers,~Pollen, summarise,Production=mean(pup.binary))
avg.probs #28% of sunflower colonies made pupae; NO buckwheat colonies made pupae!

p0<-ggplot(data=avg.probs, aes(Pollen, Production)) + 
    geom_bar(colour="black", stat="identity", size=2) +
  scale_x_discrete(breaks=levels(avg.probs$Pollen), 
                   labels=c("Buckwheat", "Sunflower"))+
  theme_classic(base_size=24)
p0 #why no axes
p001<-p0 + theme (axis.line.y=element_line(color="black", size=2),
            axis.line.x=element_line(color="black", size=2))
p01<- p001 +
  ylab("Proportion with pupae")  + xlab('Pollen') 
p01
p02<-p01 + theme(text=element_text(face="bold"))
p02

#remove space?
p03<-p02 + scale_y_continuous(expand = c(0, 0), limits = c(0, .6)) 
p03

setwd("C:/Users/Evan/Dropbox/Jess Sum15/jess stats/colony.performance")
pdf("pupae.blocks.1and2.pdf", 
    height=6, width=5, paper="special")
p03
dev.off()
#what about using bayesian framework?
###extra: Ben Bolker on adding weak priors
#http://stats.stackexchange.com/questions/132677/binomial-glmm-with-a-categorical-variable-with-full-successes
#install.packages("blme")
library(blme)
blme_ <- bglmer(pup.binary~Pollen + Infection+(1|Colony),
                data=Perf03, family=binomial,
                       fixef.prior = normal(cov = diag(9,3)))
#9 is for variance of 3 on logit scale (very large), 
#3 is for 3 terms in model (intercept, pollen, infection)

summary(blme_)


Boxplot(eggnum~Pollen*Infection, data=Perf03)
#interesting-- infection INCREASED egg production in sunflower
#but CRATERED egg production in buckwheat!!
egg1<-glmer( eggnum~ Pollen  + (1|Colony), data=Perf03,
       nAGQ = 25,
       family=poisson)
summary(egg1)
egg2<-glmer( eggnum~ Pollen + Infection  + (1|Colony), data=Perf03,
             nAGQ = 25,
             family=poisson)
summary(egg2)
egg3<-glmer( eggnum~ Pollen * Infection  + (1|Colony), data=Perf03,
             nAGQ = 25,
             family=poisson)
summary(egg3)
Anova(egg3)
library(lsmeans)
margmeans<-summary(lsmeans(object=egg3, ~Infection*Pollen))
margmeans2<-summary(lsmeans(object=egg3, ~Infection*Pollen), type="response")
margmeans2
margmeans
margmeans$btmean<-exp(margmeans$lsmean)
margmeans$rawhi<-margmeans$lsmean + margmeans$SE
margmeans$rawlo<-margmeans$lsmean - margmeans$SE
margmeans$bthi<-exp(margmeans$rawhi)
margmeans$btlo<-exp(margmeans$rawlo)
margmeans$btmean #matches margmeans2, good

