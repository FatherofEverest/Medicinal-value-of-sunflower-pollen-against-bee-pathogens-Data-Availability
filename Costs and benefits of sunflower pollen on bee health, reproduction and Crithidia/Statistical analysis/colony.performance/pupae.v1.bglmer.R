####2016.03#########
#analysis of sunflower microcolony performance

#####Pupal production################
####Source data: use summary created by LS Adler
#EPY added back colonies from which new adults emerged -- 
#colonies 3, 11, 23, 25, 33
#coded as "1" for pupal production

#rm(list=ls()) #clear memory
library(boot) #inverse logit
library(plyr)
library(dplyr)
library(car)
library(nlme)
library(lme4)
library(lsmeans)
library(glmmTMB) #takes some work to install this
library(bbmle)
library(ggplot2)
library(brglm)
library(blme)

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


##With random effect
GLMM.12<- glmer(pup.binary ~ Pollen*Infection +  (1 | Colony), data = Producers,
      family = binomial(link = "logit"), nAGQ = 100)  # Set nAGQ to # of desired iterations
Anova(GLMM.12)

#Go to bayesian model to address complete separation
#what about using bayesian framework?
###extra: Ben Bolker on adding weak priors
#http://stats.stackexchange.com/questions/132677/binomial-glmm-with-a-categorical-variable-with-full-successes

bmod1 <- bglmer(pup.binary ~ Pollen*Infection +  (1 | Colony), data = Producers, family = binomial,
                fixef.prior=normal)
Anova(bmod1)
# Chisq Df Pr(>Chisq)   
# Pollen           7.5891  1   0.005872 **
#   Infection        0.1514  1   0.697161   
# Pollen:Infection 0.0964  1   0.756193  
drop1(bmod1, test = "Chisq")

#Remove interaction
bmod2<-update(bmod1, ~. - Pollen:Infection)
drop1(bmod2, test = "Chisq")
# pup.binary ~ Pollen + Infection + (1 | Colony)
# Df    AIC     LRT   Pr(Chi)    
# <none>       53.660                      
# Pollen     1 67.661 16.0008 6.331e-05 ***
#   Infection  1 51.817  0.1573    0.6916   

Anova(bmod2)
# Chisq Df Pr(>Chisq)   
# Pollen    7.8914  1   0.004967 **
#   Infection 0.1547  1   0.694124 
#same as for drop1()


################### Plots #############
Lsm<-lsmeans(bmod2, ~Pollen)
Pairs<-contrast(Lsm, "pairwise") 
Pairs
# contrast  estimate       SE df z.ratio p.value
# B - S    -3.403488 1.211569 NA  -2.809  0.0050

##Sunflowers distinguish themselves

##Get letters to plot
Lets<-cld(Lsm, Letters = letters)
Lets
Lets$.group<-gsub(" ", "", Lets$.group) ##remove spaces

################## PLOT ON SCALE OF RESPONSE ##########
#Let's use inverse logit in boot to get original mean+SE and mean-SE
Lsmeans.df<-as.data.frame(summary(lsmeans(bmod2, ~Pollen)))
Lsmeans.df$pup.binary<-inv.logit(Lsmeans.df$lsmean)
Lsmeans.df$btlo<-inv.logit(Lsmeans.df$lsmean - Lsmeans.df$SE) #lower bound
Lsmeans.df$bthi<-inv.logit(Lsmeans.df$lsmean + Lsmeans.df$SE) #upper bound

Lsmeans.df

#Merge with letters
Lets2<- dplyr::select(Lets, Pollen, .group)
Lets2

Lsmeans.df.merged<- join(Lsmeans.df, Lets2, by = "Pollen")

Lsmeans.df.merged
# Pollen    lsmean       SE df asymp.LCL  asymp.UCL pup.binary       btlo       bthi .group
# 1      B -5.058268 1.739015 NA -8.466674 -1.6498613 0.00631641 0.00111556 0.03491657      a
# 2      S -1.654780 1.328725 NA -4.259034  0.9494738 0.16046395 0.04817663 0.41920088      b
#3x as many larvae

Lsmeans.df<-Lsmeans.df.merged #overwrites old

#x-axis:
pollen.axis<-c("Buck", "Sun")
#y-axis:
ylabel<-"Proportion with pupae" 

#ready to plot?
p<-ggplot(Lsmeans.df, aes(x = Pollen, y = pup.binary)) + 
  # #Option to sprinkle raw data... 
  # geom_point(data = Perf03, aes( x = Pollen, y = pup.binary), position = "jitter", 
  #            color = "blue", shape = 25, size = 0.5, alpha = 0.4)+
  geom_pointrange(aes(ymin = btlo, ymax = bthi), size = 1, color = "black") +
  scale_x_discrete(labels=pollen.axis) +
  ylab(ylabel) + #y label
  xlab("Pollen")
p
#Add text
Pfig<- p + geom_text(aes(y = bthi + 0.2 * max(bthi), label = .group),
                     size = 8)
Pfig

##Export
# setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/jess stats/colony.performance")
# ggsave("pupae.blocks.1and2.v2.bglmer.pdf", height = 5, width =5)

##Consolidate metrics
#A. Number of larvae
#("Larvalnumber.v3.TMB.R")
Lfig
#B. Mean larval mass
#("Larvalmass.v1.lme4.R")
P.lmass
#C. Proportion with pupae
Pfig

#Arrange and plot
plot_grid(Lfig, P.lmass, Pfig, nrow=1, labels = LETTERS,
          label_size = rel(20))
ggsave("MC.performance.3panel.pdf", height = 4, width = 10)

###extra: Ben Bolker on adding weak priors
#http://stats.stackexchange.com/questions/132677/binomial-glmm-with-a-categorical-variable-with-full-successes
#install.packages("blme")
library(blme)
blme_ <- bglmer(pup.binary~Pollen + Infection+(1|Colony),
                data=Perf03, family=binomial,
                       fixef.prior = normal(cov = diag(9,3)))
#9 is for variance of 3 on logit scale (very large), 
#3 is for 3 terms in model (intercept, pollen, infection)


