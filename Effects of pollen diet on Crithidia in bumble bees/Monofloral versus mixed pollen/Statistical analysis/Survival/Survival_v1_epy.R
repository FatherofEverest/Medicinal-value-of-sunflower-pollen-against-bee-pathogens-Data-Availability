#Mixed model survival analysis
#Individual bees fed different pollens over 7d
library(coxme)
library(survival) 
library(car)
library(ggplot2)
#google drive (epy)
setwd("C:/Users/Evan/Google Drive/ The Sunflower Pollen Project/Analysis.Individual.Bee.and.Farm/Analyses/Survival")
survivedata<-read.csv('Pollen Exp. Data_Mortality.csv')

#jon's machine:
survivedata<-read.csv("/Users/rg9403/Desktop/PollenEXP/PollenExp_Winter2014/DATA/CSV/Pollen Exp. Data_Mortality.csv", header=TRUE)

str(survivedata)
View(survivedata)
head(survivedata)
attach(survivedata)

library(plyr)
library(dplyr)
aggregate(survivedata$Dead.Before.Dissect~survivedata$Treatment, FUN=mean)


# let's set up the variables
time<-Time.To.Death 
event<-Dead.Before.Dissect
X<-Treatment
time<-as.numeric(paste(time))
X<-as.factor(X)
library(coxme)


# Let's check the structure of the data, everything looks good
str(survivedata$Treatment)
str(survivedata)

str(survivedata$Time.To.Death)
#int [1:250] 7 3 7 7 7 7 7 1 0 1 ...


#### Cox-mixed effects model ####
MIXED0 <- coxme(Surv(Time.To.Death,Dead.Before.Dissect) ~ Treatment + Callow.Mass + (1|Inoc.Date) + (1|Colony.ID), data=survivedata)
summary(MIXED0) #no error, no treat effects
library(car)
Anova(MIXED0)

#drop mass
M1<-update(MIXED0, ~. - Callow.Mass)
summary(M1)
Anova(M1)

#let's see if the random effects are helpful
FULL<-coxme(Surv(Time.To.Death,Dead.Before.Dissect) ~ Treatment  + (1|Inoc.Date) + (1|Colony.ID), data=survivedata)
NOCOL<- coxme(Surv(Time.To.Death,Dead.Before.Dissect) ~ Treatment  + (1|Inoc.Date) , data=survivedata)
NODATE<-coxme(Surv(Time.To.Death,Dead.Before.Dissect) ~ Treatment  + (1|Colony.ID), data=survivedata)
Fixed<-coxph(Surv(Time.To.Death,Dead.Before.Dissect) ~ Treatment  , data=survivedata)

library(bbmle)
AICtab(FULL, NOCOL, NODATE, Fixed)
#inoculation date could be dropped

anova(FULL, NOCOL) #keep colony random effect
anova(FULL, NODATE) #drop inoculation date random effect

#NODATE is our reduced model
summary(NODATE)
Anova(NODATE) #treatment not significant (chi square test)
Stripped<-update(NODATE, ~. - Treatment)
anova(NODATE, Stripped) #treatment not significant, Likelihood ratio test


kmsurvival1 <- survfit(Surv(time, event) ~ X)
summary(kmsurvival1)
# Call: survfit(formula = Surv(time, event) ~ X)
# 
# X=Buck 
# time n.risk n.event survival std.err lower 95% CI upper 95% CI
# 0     62       2    0.968  0.0224        0.925        1.000
# 2     60       1    0.952  0.0273        0.900        1.000
# 3     59       1    0.935  0.0312        0.876        0.999
# 6     58       4    0.871  0.0426        0.791        0.959
# 
# X=Mix 
# time n.risk n.event survival std.err lower 95% CI upper 95% CI
# 0     60       1    0.983  0.0165        0.951        1.000
# 1     59       1    0.967  0.0232        0.922        1.000
# 3     58       3    0.917  0.0357        0.849        0.989
# 4     55       2    0.883  0.0414        0.806        0.968
# 5     53       2    0.850  0.0461        0.764        0.945
# 6     51       1    0.833  0.0481        0.744        0.933
# 7     50       1    0.817  0.0500        0.724        0.921
# 
# X=Rape 
# time n.risk n.event survival std.err lower 95% CI upper 95% CI
# 1     62       1    0.984  0.0160        0.953        1.000
# 2     61       2    0.952  0.0273        0.900        1.000
# 3     59       1    0.935  0.0312        0.876        0.999
# 5     58       5    0.855  0.0447        0.772        0.947
# 6     53       3    0.806  0.0502        0.714        0.911
# 7     50       1    0.790  0.0517        0.695        0.898
# 
# X=Sun 
# time n.risk n.event survival std.err lower 95% CI upper 95% CI
# 0     66       2    0.970  0.0211        0.929        1.000
# 1     64       2    0.939  0.0294        0.884        0.999
# 5     62       2    0.909  0.0354        0.842        0.981

plot(kmsurvival1, xlab="Time (d)", ylab="Survival Probability")
levels(X)
#[1] "Buck" "Mix"  "Rape" "Sun"


#jpeg("survival_v1HOR.jpg") #uncomment for jpeg
#pdf("survival_v1HOR.pdf", height=7, width=6, paper="special") #uncomment for pdf
par(mfrow=c(1,1))
par(mar=c(5,7,1,1)) #bottom, left, top right-- note biggest margin on left 
plot(kmsurvival1, ylim=c(0.75, 1.02), xlab="Time (d)", 
     ylab="Proportion Surviving\n", #added 
     cex.lab=2, lty=c(1,5,3,4), mark.time=c(7), mark=c(0,21,24,25), cex=2,
     frame.plot=FALSE, #no border
     cex.axis=1.5, las=1) #las=1 : all text horizontal
legend(0.5, 0.9, levels(X), cex=2, lty=c(1,5,3,4), pch=c(1,21,24,25), 
       title = "Pollen")
dev.off()


