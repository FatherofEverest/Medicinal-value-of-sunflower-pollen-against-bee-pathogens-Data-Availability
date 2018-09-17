library(survival) 

survivedata<-read.csv("/Users/rg9403/Desktop/PollenEXP/PollenExp_Winter2014/DATA/CSV/Pollen Exp. Data_Mortality.csv", header=TRUE)

View(survivedata)
head(survivedata)
attach(survivedata)

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
# 'data.frame':	250 obs. of  10 variables:
#   $ Bee........ID      : int  1 3 4 5 6 7 8 9 10 11 ...
# $ Colony.ID          : Factor w/ 6 levels "JG 1","JG 3",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Treatment          : Factor w/ 4 levels "Buck","Mix","Rape",..: 4 3 2 4 1 3 2 4 1 3 ...
# $ Date.Emerged       : Factor w/ 39 levels "01-Feb","01-Jan",..: 36 36 36 36 38 38 38 38 38 38 ...
# $ Callow.Mass        : num  0.15 0.101 0.151 0.097 0.161 ...
# $ Inoc.Date          : Factor w/ 37 levels "01-Feb","01-Jan",..: 2 2 2 2 4 4 4 4 4 4 ...
# $ Inoc..Day.of.Week  : Factor w/ 8 levels "Fri","Mon","Sat",..: 6 6 6 6 1 1 1 1 1 1 ...
# $ DOD.Weight..g.     : Factor w/ 197 levels ".","0.0577","0.0602",..: 131 1 70 26 106 42 18 1 1 1 ...
# $ Dead.Before.Dissect: int  0 1 0 0 0 0 0 1 1 1 ...
# $ Time.To.Death      : int  7 3 7 7 7 7 7 1 0 1 ...
str(survivedata$Time.To.Death)
#int [1:250] 7 3 7 7 7 7 7 1 0 1 ...


#### Cox-mixed effects model ####
coxme <- coxme(Surv(Time.To.Death,Dead.Before.Dissect) ~ Treatment + Callow.Mass + (1|Inoc.Date) + (1|Colony.ID), data=survivedata)
summary(coxme) #no error, no treat effects
# Cox mixed-effects model fit by maximum likelihood
# Data: survivedata
# events, n = 38, 250
# Iterations= 6 34 
# NULL Integrated   Fitted
# Log-likelihood -206.8513  -199.9669 -196.078
# 
# Chisq   df        p  AIC   BIC
# Integrated loglik 13.77 6.00 0.032328 1.77 -8.06
# Penalized loglik 21.55 7.04 0.003134 7.46 -4.08
# 
# Model:  Surv(Time.To.Death, Dead.Before.Dissect) ~ Treatment + Callow.Mass +      (1 | Inoc.Date) + (1 | Colony.ID) 
# Fixed coefficients
# coef   exp(coef)  se(coef)     z    p
# TreatmentMix   0.3782626 1.459746184 0.4650516  0.81 0.42
# TreatmentRape  0.4574683 1.580068688 0.4521777  1.01 0.31
# TreatmentSun  -0.4530057 0.635714494 0.5434546 -0.83 0.40
# Callow.Mass   -5.2264640 0.005372489 4.2445770 -1.23 0.22
# 
# Random effects
# Group     Variable  Std Dev      Variance    
# Inoc.Date Intercept 0.0199483113 0.0003979351
# Colony.ID Intercept 0.5975996711 0.3571253669



library(car)
#Anova in "car" only takes certain kinds of model object
Anova(coxme) # Didn't work
#you should be able to compare a model that does not have treatment
#with the full model
coxme.notreat<-coxme(Surv(Time.To.Death,Dead.Before.Dissect) ~  Callow.Mass + (1|Inoc.Date) + (1|Colony.ID), data=survivedata)
#now anova command to compare full model vs model with no pollen treat
#anova(fullmodel, reducedmodel)
anova(coxme, coxme.notreat)
#"to assess significance of treatment effect on mortality, we compared log-likelihood of #nested models with and without diet treatment as a predictor using a chi-squared test"
# Analysis of Deviance Table
# Cox model: response is  Surv(Time.To.Death, Dead.Before.Dissect)
# Model 1: ~Treatment + Callow.Mass + (1 | Inoc.Date) + (1 | Colony.ID)
# Model 2: ~Callow.Mass + (1 | Inoc.Date) + (1 | Colony.ID)
# loglik  Chisq Df P(>|Chi|)
# 1 -199.97                    
# 2 -202.21 4.4906  3    0.2131
 #p=0.21 for effect of treatment
#chi-squared test: chisq 4.49  on 3 df (with 4 treatments)
#you can stop there -- 



#alternative is to do a wald test where you remove the fixed coefficient, 
#see the script using aod, do wald.test (model, terms=c(1:3))
#this is chi squared test for individual predictors or groups of them
#ie your different diet treats
#in theory these are both chisquare tests
#but they might be comparing different things-- ie wald test different than log-likelihood based aov comparison
library(aod)
wald.test(vcov(coxme),fixef(coxme), Terms=c(1:3))
# Wald test:
#   ----------
#   
#   Chi-squared test:
#   X2 = 4.1, df = 3, P(> X2) = 0.25
#result is close, i think it's a little different because
#one of the diet treatments is reference level,
#so we only took away 3 terms instead of 4
#if you use this:
########"used a Wald test to examine the significance of diet treatment in the mixed model"

#you can also try to simplify your model a bit
#using the same algorithm as above but test mass
wald.test(vcov(coxme),fixef(coxme), Terms=c(4)) #mass p=0.22
#you could make a new coxme w/o mass, if desired, 
#then do the same nested model aov w/ and w/o treat
#and also the wald test removing terms 1,2,3



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
par(mfrow=c(1,1))
plot(kmsurvival1, ylim=c(0.7, 1.02), xlab="Time (d)", ylab="Survival Probability", cex.lab=1.5, lty=c(1,5,3,4), mark.time=c(7), mark=c(0,21,24,25), cex=1.6)
legend(0.6, 0.87, levels(X), cex=0.95, lty=c(1,5,3,4), pch=c(1,21,24,25), title = "pollen")

