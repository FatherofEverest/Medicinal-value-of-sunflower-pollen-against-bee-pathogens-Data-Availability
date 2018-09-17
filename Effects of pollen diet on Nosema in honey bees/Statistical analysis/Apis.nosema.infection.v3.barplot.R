############################################################################################################################################################
# Honey Bee - Nosema Exp Analysis
# Performed by JJG - May 25, 2017
 ##Edited EPY May 27 to add random effect of cup and fixed effect of time
#updated 10 June 2017 to use bar chart

#Groups (BeeCupID) of 50 Honey bees (Apis mellifera) were inoculated with Nosema spores and provided either
#sunflower pollen, buckwheat pollen, or no pollen (as a control diet).  Groups of 5 bees and 10 bees were sacraficed on days
#10 and 15, repectively, to measure Nosema infection intensity. 

#Model Spores ~ Treatment * Time + (1|Cup)

#Data wrangling
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)

#Exploratory and "Anova"
library(MASS)
library(car)

#Models
library(glmmADMB)
library(glmmTMB)
library(lme4)
library(lsmeans)

#Model Extensions & graphs
library(lsmeans)
library(ggplot2)
library(cowplot)


#EPY
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Analysis.Individual.Bee.and.Farm/nosema.apis")
data1<-read.csv("Honey Bee Nosema Exp_June 2016 (2).csv", header=TRUE)
source("glmmTMB.to.lsmeans.R") #helper function for TMB to talk to lsmeans

#NCSU:
#data1<-read.csv("/Users/rg9403/Google Drive/Sunflower Pollen Experiments/Irwin Lab_Sunflower Pollen Experiments/Exp 6 _ Apis mellifera_Nosema_Sun_Buck_NoPollen_June 2016/Honey Bee Nosema Exp_June 2016 (2).csv", header=TRUE)

str(data1)

#Rename:
D<-dplyr::rename(data1, Cup = BeeCupID, 
          Count.10d = NosemaDay10.Count,
          Count.15d =Average.Day.15.Count)
D<-dplyr::select(D, Cup, Treatment, Count.10d, Count.15d)
##Melt to long format
Dm<-melt(data = D, id.vars = c("Cup", "Treatment"), 
         measure.vars = c("Count.10d", "Count.15d"),
         variable.name = "Time",
         value.name = "Count")
str(Dm)
str(Dm$Time)
Dm$Time<-plyr::mapvalues(Dm$Time, 
                         from = c("Count.10d", "Count.15d"),
                         to = c("D10" ,"D15"))

###Response variable: Exploratory (Distributions)

densityPlot(Dm$Count) #bimodal
car::Boxplot(Dm$Count ~ Dm$Treatment + Dm$Time )

qqp(Dm$Count, "norm") #super
### Not too shabby
poisson <- fitdistr(Dm$Count, "Poisson")
poisson
poisson$estimate
?qqp
qqp(x = Dm$Count, distribution = "pois", lambda = poisson$estimate)
# Looks nearly the same as the Normal fit

nbinomDist <- fitdistr(Dm$Count, "Negative Binomial")
qqp(x = Dm$Count, distribution = "nbinom",
    size = nbinomDist$estimate[[1]],mu = nbinomDist$estimate[[2]])

### ERRORS: Error in solve.default(res$hessian) : 
###Lapack routine dgesv: system is exactly singular: U[2,2] = 0
  ##EPY: no errors for me

##############################################################################
# Nosema Infection intensity models - GLMs
   #EPY edit: GLMM with Cup as random effect
##############################################################################

#Relevel treatments: 
levels(Dm$Treatment)
Dm$Treatment<-factor(Dm$Treatment, levels = c("NP", "B", "S"))

#Use Poisson distribution 

FitPoisson<- glmmTMB(Count~Treatment*Time + (1|Cup), 
                     family = poisson(), data = Dm)
#Much faster than ADMB!!

summary(FitPoisson) #overdispersion because residual deviance much larger than df
# AIC      BIC   logLik deviance df.resid 
# 2176.5   2192.1  -1081.3   2162.5       61 

#Compare to negative binomial
FitNB<-glmmTMB(Count~Treatment*Time + (1|Cup), 
               family = "nbinom2", data = Dm)
summary(FitNB)
# AIC      BIC   logLik deviance df.resid 
# 946.0    963.7   -465.0    930.0       60 

# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          7.42953    0.06161  120.59  < 2e-16 ***
#   TreatmentNP         -0.72043    0.08562   -8.41  < 2e-16 ***
#   TreatmentS          -0.41616    0.08806   -4.73 2.29e-06 ***
#   TimeD15              0.19262    0.07851    2.45   0.0142 *  
#   TreatmentNP:TimeD15 -0.78611    0.10973   -7.16 7.84e-13 ***
#   TreatmentS:TimeD15   0.14954    0.11154    1.34   0.1800    


library(bbmle)
AICtab(FitPoisson, FitNB)
# dAIC   df
# FitNB         0.0 8 
# FitPoisson 1230.6 7 
#HUGE change in AIC, it seems

# Let's see if any of the models better than the others
anova(FitPoisson, FitNB, test="Chisq") 
#Df     AIC     BIC   logLik deviance  Chisq Chi Df Pr(>Chisq)    
# FitPoisson  7 2176.53 2192.07 -1081.27  2162.53                             
# FitNB       8  945.97  963.72  -464.98   929.97 1232.6      1  < 2.2e-16 ***
 #Stick with FitNB

summary(FitNB)
Anova(FitNB) #unhappy, Anova not updated to deal with glmmTMB?
drop1(FitNB, test = "Chisq")
# Df    AIC    LRT  Pr(>Chi)    
# <none>            945.97                     
# Treatment:Time  2 990.19 48.226 3.371e-11 ***



plot(residuals(FitNB)) #that's from heaven

#Try ADMB fit instead
Dm$Cup<-as.factor(Dm$Cup)
FitNB.A<-glmmadmb(Count~Treatment*Time + (1|Cup), 
               family = "nbinom", data = Dm)
summary(FitNB.A)
Anova(FitNB.A)
# Response: Count
# Df   Chisq Pr(>Chisq)    
# Treatment       2 92.3691  < 2.2e-16 ***
#   Time            1  6.3682    0.01162 *  
#   Treatment:Time  2 52.9712  3.144e-12 ***

#compare glmmadmb against glmmTMB
summary(FitNB)
summary(FitNB.A) #perfect match!
  #Looks like we have to drop1() for significance tests and Ben Bolker helper function for lsmeans
   #from glmmTMB
    #but the model coefficients are the same

##Pairwise comparisons
library(lsmeans)
Lsm.nb<-lsmeans(FitNB, ~Treatment|Time)
nb.contrasts<- contrast(Lsm.nb, method ="pairwise")
nb.contrasts
# Time = D10:
#   contrast   estimate         SE df z.ratio p.value
# B - NP    0.7204282 0.08561876 NA   8.414  <.0001
# B - S     0.4161561 0.08805502 NA   4.726  <.0001
# NP - S   -0.3042721 0.08656078 NA  -3.515  0.0013
# 
# Time = D15:
#   contrast   estimate         SE df z.ratio p.value
# B - NP    1.5065399 0.08611563 NA  17.494  <.0001
# B - S     0.2666119 0.08710224 NA   3.061  0.0062
# NP - S   -1.2399280 0.08614269 NA -14.394  <.0001

#No Pollen lowest, then sunflower, then buckwheat
#Compare the two timepoints
Lsm.treatment<- lsmeans(FitNB, ~Treatment)
contrast(Lsm.treatment, "pairwise")
# contrast   estimate         SE df z.ratio p.value
# NP - B   -1.1134843 0.06605224 NA -16.858  <.0001
# NP - S   -0.7721019 0.06636456 NA -11.634  <.0001
# B - S     0.3413824 0.06752839 NA   5.055  <.0001
Lsm.treatment.df<-as.data.frame(summary(Lsm.treatment))
Lsm.treatment.df
Lsm.treatment.df$Count<-exp(Lsm.treatment.df$lsmean)
Lsm.treatment.df
# Treatment   lsmean         SE df asymp.LCL asymp.UCL     Count
# 1        NP 6.412356 0.04597292 NA  6.322251  6.502462  609.3277
# 2         B 7.525841 0.04741854 NA  7.432902  7.618779 1855.3721
# 3         S 7.184458 0.04801578 NA  7.090349  7.278567 1318.7745

#Reduction in sunflower compared to buckwheat
Reduction.sun<-(Lsm.treatment.df$Count[Lsm.treatment.df$Treatment=="B"]-
              Lsm.treatment.df$Count[Lsm.treatment.df$Treatment=="S"])/
  Lsm.treatment.df$Count[Lsm.treatment.df$Treatment=="B"]

Reduction.sun
# > Reduction.sun
# [1] 0.2892129

#Sunflower vs no pollen
Reduction.np.vs.sun<-(Lsm.treatment.df$Count[Lsm.treatment.df$Treatment=="S"]-
                        Lsm.treatment.df$Count[Lsm.treatment.df$Treatment=="NP"])/
  Lsm.treatment.df$Count[Lsm.treatment.df$Treatment=="S"]

Reduction.np.vs.sun
#[1] 0.5379591

#No pollen vs buckwheat
Reduction.np.vs.b<-(Lsm.treatment.df$Count[Lsm.treatment.df$Treatment=="B"]-
                        Lsm.treatment.df$Count[Lsm.treatment.df$Treatment=="NP"])/
  Lsm.treatment.df$Count[Lsm.treatment.df$Treatment=="B"]
Reduction.np.vs.b
#0.6715873

#Analysis by time
Lsm.time<- lsmeans(FitNB, ~Time)
contrast(Lsm.time, "pairwise")
# contrast    estimate         SE df z.ratio p.value
# D10 - D15 0.01957085 0.04512284 NA   0.434  0.6645
Lsm.time.df<-as.data.frame(summary(Lsm.time))
Lsm.time.df
Lsm.time.df$Count<-exp(Lsm.time.df$lsmean)
Lsm.time.df
# Time  lsmean         SE df asymp.LCL asymp.UCL    Count
# 1  D10 7.05067 0.03543993 NA  6.981209  7.120131 1153.632
# 2  D15 7.03110 0.03529116 NA  6.961930  7.100269 1131.274


Lsmeans.df<- as.data.frame(summary(Lsm.nb))
#significance letters
Lets<-cld(Lsm.nb, method = "pairwise", Letters = letters)
Lets
Lets<-as.data.frame(Lets)
Lets
Lets.slim<-dplyr::select(Lets, Treatment, Time, .group)
Lsm.w.lets<-join(Lets.slim, Lsmeans.df, by = c("Treatment", "Time"))
Lsm.w.lets

Lsm.w.lets$.group <-
  gsub(pattern= " ", replacement = "", x= Lsm.w.lets$.group, 
       ignore.case = FALSE, perl = FALSE,
     fixed = FALSE, useBytes = FALSE)
Lsm.w.lets

#Option: add a space before each letter to nudge to right on graph
#Lsm.w.lets$.group<- paste(" ", Lsm.w.lets$.group)

Lsmeans.df<-Lsm.w.lets

Lsmeans.df$btmean<-exp(Lsmeans.df$lsmean)
Lsmeans.df$btlo<-exp(Lsmeans.df$lsmean - Lsmeans.df$SE) #lower bound
Lsmeans.df$bthi<-exp(Lsmeans.df$lsmean + Lsmeans.df$SE) #upper bound

Lsmeans.df
# Treatment Time .group   lsmean         SE df asymp.LCL asymp.UCL    btmean      btlo      bthi
# 1        NP  D10      a 6.709104 0.05945879 NA  6.592567  6.825641  819.8355  772.5100  870.0603
# 2         S  D10      b 7.013377 0.06294891 NA  6.889999  7.136754 1111.4012 1043.5962 1183.6116
# 3         B  D10      c 7.429531 0.06160947 NA  7.308778  7.550283 1685.0167 1584.3370 1792.0943
# 4        NP  D15      a 6.115609 0.06022761 NA  5.997565  6.233653  452.8717  426.4014  480.9852
# 5         S  D15      b 7.355540 0.06164379 NA  7.234720  7.476359 1564.8411 1471.2913 1664.3390
# 6         B  D15      c 7.622150 0.06151016 NA  7.501593  7.742708 2042.9504 1921.0750 2172.5579

#Labels and treatment names
ylabel<-expression(italic(Nosema)~count~("cells*"~0.02~mu*L^-1)) #"*" for no-space
xlabel<- "Time (d)"
x.ticklabs<- c("10", "15")
Pollen.names<- c("No Pollen", "Buckwheat", "Sunflower")
pd<-position_dodge(width=0.9)
plot0<- ggplot(Lsmeans.df, aes(x=Time, y=btmean, group = Treatment)) + 
  geom_bar(aes(fill=Treatment), position=pd, stat="identity", color = "black") +
  geom_errorbar(aes(ymin=btlo, ymax=bthi),
                size=1, width=0.4,  # Width of the error bars
                position=pd) +
  theme_cowplot() 
plot0

p1<- plot0 + ylab (ylabel) + xlab (xlabel) +
  scale_x_discrete(labels = x.ticklabs) + 
  scale_fill_manual(name = NULL, labels = Pollen.names, 
                    values = c("white", "blue", "orange"))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.2 * max(Lsmeans.df$bthi)))+
  theme(axis.line = element_line(size = 2),
        text = element_text (face = "bold"))
p1

#Bonus 1: Add significance letters
p.let<- p1 + geom_text(aes(y = bthi + 0.1* max(bthi), label = .group), 
                       position=pd, 
                       color = "black",
                       size = 8)
p.let

getwd()
# ggsave("Nosema.sunflower.v2.barchart.pdf", height = 4, width = 6)

#Proceed to infection script to make panel B
#Nosema Exp Survival Analysis.v2.ggsurv.R

