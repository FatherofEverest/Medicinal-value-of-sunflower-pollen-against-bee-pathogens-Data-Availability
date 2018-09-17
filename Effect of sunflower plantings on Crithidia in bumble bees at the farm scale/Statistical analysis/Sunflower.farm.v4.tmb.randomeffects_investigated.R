#######SCRIPT TO ANALYZE CRITHIDIA FROM WILD BEES########
###Collected from Amherst area farms summer 2015
####PROPOSED EXPLANATORY VARIABLE: sunflower area#

#Updated 2017.05.26 -- lsmeans behaves strangely; feed glmmTMB model instead

#Updated 2018.02.06:
##For R1: explicit test of Farm
# L1.nofarm <- update(L1plus.tmb, ~. - (1|Farm))
# anova(L1.nofarm, L1plus.tmb)
##"Farm" was not a meaningful explainer of variance in infection intensity


rm(list=ls())
#install.packages("glmmADMB", repos="http://r-forge.r-project.org",type="source")
library(glmmADMB)
library(glmmTMB)
library(lsmeans)
library(ggplot2)
library(emmeans)##new version of lsmeans

#read data
#Google Drive (epy)
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Analysis.Individual.Bee.and.Farm/Analyses/Farm Sampling")
farmdata<-read.csv("Analysis_1_All.csv")

#JG machine:
farmdata<-read.csv("/Users/rg9403/Desktop/sunflowersamplepatrickdata/SUNFLOWER_FARM_SAMPLING/Sunflower Sampling Data_CSVs/Analysis_1_All.csv",header=TRUE)
str(farmdata)
# 'data.frame':	667 obs. of  13 variables:
#   $ CollectionDate : int  : Julian date of collection
# $ Farm           : Factor w/ 22 levels : denotes farm where bees collected
# $ BeeID  : unique number for each bee      
# $ WingSize       : num  : size of marginal cell of right forewing
# $ FlowerID       : Factor w/ 18 levels "BAS","BRG","BSUN",..: 2 3 3 3 3 3 7 9 9 9 ...
# $ Count          : int  : count of motile Crithidia cells in 0.02uL of gut supernatant
# $ Log1P.Count.   : log10 of "Count"?
# $ FlowerID2 ????means what???     : Factor w/ 2 levels "NONSUN","SUN": 1 2 2 2 2 2 1 1 1 1 ...
# $ FlowerID3 ???means what?     : Factor w/ 3 levels ".","HELI","NONHELI": 3 2 2 2 2 2 3 3 3 3 ...
# $ FarmSize  ???units??     : Factor w/ 22 levels "112,137.41","131,041.50",..: 10 10 10 10 10 10 10 10 10 10 ...
# $ Log1P.FarmSize. ???units???: num  4.52 4.52 4.52 4.52 4.52 4.52 4.52 4.52 4.52 4.52 ...
# $ SunArea        : num  67.9 67.9 67.9 67.9 67.9 ...
# $ Log1P.SunArea. : 
###     EPY COMMENT-- I HAVE A PROBLEM WITH "Log1P.SunArea"
###There are zeroes for "SunArea"

farmdata$Count<-as.numeric(farmdata$Count)

Data<-na.omit(farmdata) #in order for glmmadmb's to run

Data1<-Data



L1plus.tmb<-glmmTMB(Count ~ Log1P.SunArea. + CollectionDate +(1|Farm), 
                 family = "nbinom2",ziformula = ~0, data=farmdata) 
#check it matches admb model
summary(L1plus.tmb)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    11.89645    3.04379   3.908 9.29e-05 ***
#   Log1P.SunArea. -0.26368    0.09956  -2.648  0.00809 ** 
#   CollectionDate -0.03758    0.01329  -2.827  0.00470 ** 
summary(L1plus) #fixed effects are same; slight difference in se(intercept)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    11.89628    2.12500    5.60  2.2e-08 ***
#   Log1P.SunArea. -0.26369    0.09749   -2.70   0.0068 ** 
#   CollectionDate -0.03758    0.00927   -4.05  5.1e-05 ***

##Get some concrete meanings of the parameters:
##How much more sunflower to change infection intensity by factor of 2?
#This means that infection on log scale must change by ln(2)
dY <- log(2) #0.6931472
str(summary(L1plus.tmb))
summary(L1plus.tmb)$coefficients$cond[1,]
Suncoef <- summary(L1plus.tmb)$coefficients$cond[2,1]
Daycoef <- summary(L1plus.tmb)$coefficients$cond[3,1]

##how much change in Crithidia per one-unit (10x) change in sunflower area?
(dCrith.per.dSun <- 1 - exp(Suncoef))
#[1] 0.2317803
#CI's:
Coef.df <- as.data.frame(summary(L1plus.tmb)$coefficients$cond)
Coef.df
Coef.df$Est.upper <- Coef.df$Estimate + 2*Coef.df$`Std. Error`
Coef.df$Est.lower <- Coef.df$Estimate - 2*Coef.df$`Std. Error`
Coef.df$Term <- row.names(Coef.df)

Df.to.exp <- dplyr::select(Coef.df, Estimate, Est.lower, Est.upper)
Df.exp <- exp(as.matrix(Df.to.exp))
Df.exp

Df.exp.delta.percent <- 100*(1 - Df.exp)
Df.exp.delta.percent
# Estimate     Est.lower     Est.upper
# (Intercept)    -1.467433e+07 -33224.217608 -6.461935e+09
# Log1P.SunArea.  2.317803e+01     37.048287  6.251705e+00
# CollectionDate  3.688502e+00      6.215458  1.093458e+00
##6.25 to 37% decrease in Crithidia per 10x increase in sunflower area

#how much d(x) needed to change absolute count by factor of 2?
#y1 = ln(z)
#y2 = ln(2z)
#dy = ln(2z) - ln(z) = ln(2z/z) = ln(2)

#Our model formula derivative is:
#dy = b * d(x)
#substitute in for dy
#dy = ln(2)
#ln(2) = b*dx
#ln(2)/b = dx

##For sunflower
(dSun <- -dY / Suncoef)
#[1] 2.628749  

(dSun.raw = 10^dSun)
#425.3527

#Our model predicts a halving of infection for every 425x increase in sunflower area


##For days: interpret the coefficient
(dDay <- -dY / Daycoef)
#18.44336 

#Our model predicts a halving of infection 
  #for every 18.4 d sampling


##For R1: explicit test of Farm
L1.nofarm <- update(L1plus.tmb, ~. - (1|Farm))
anova(L1.nofarm, L1plus.tmb)
##"Farm" was not a meaningful explainer of variance in infection intensity
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# L1.nofarm   4 3893.7 3911.7 -1942.9   3885.7                         
# L1plus.tmb  5 3895.5 3918.0 -1942.8   3885.5 0.2406      1     0.6238

summary(L1.nofarm)
summary(L1plus.tmb)

# The random effect explained how much variance?
#   Let's compare residuals of Model without farm to model with farm as random effect
str(L1plus.tmb)
##Coefficients for Farm
# Groups Name        Variance Std.Dev.
# Farm   (Intercept) 0.03346  0.1829  

#Try Bolker suggestion to estimate r-squared -- compare residuals of model vs intercept-only model
#https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#storing-information
#1-var(residuals(m))/(var(model.response(model.frame(m)))
m <- L1plus.tmb
rsq.full <- 1-var(residuals(m))/(var(model.response(model.frame(m))))
#0.1301256
m2 <- L1.nofarm
rsq.reduced <- 1-var(residuals(m2))/(var(model.response(model.frame(m2))))
#0.1134202
rsq.full - rsq.reduced
#0.01670544

#Get stats for fixed effects
drop1(L1plus.tmb, test = "Chisq")
# Count ~ Log1P.SunArea. + CollectionDate + (1 | Farm)
# Df    AIC     LRT  Pr(>Chi)    
# <none>            3895.5                      
# Log1P.SunArea.  1 3900.4  6.8794 0.0087195 ** 
#   CollectionDate  1 3908.0 14.4663 0.0001427 ***
MyModel<-L1plus.tmb


areas<-seq(-0.2,max(0.1+farmdata$Log1P.SunArea.), 0.1)

#script to allow lsmeans to talk with glmmTMB
source("glmmTMB.to.lsmeans.R")



########marginal means on scale of linear predictor##########
margmeans<-summary(lsmeans(MyModel, ~Log1P.SunArea., 
                           at=list(Log1P.SunArea.=areas)))
#Manually exponentiate the (mean+-SE) from scale of linear predictor
margmeans$btmean<-exp(margmeans$lsmean)
margmeans$rawhi<-margmeans$lsmean + margmeans$SE
margmeans$rawlo<-margmeans$lsmean - margmeans$SE
margmeans$bthi<-exp(margmeans$rawhi)
margmeans$btlo<-exp(margmeans$rawlo)                        
margmeans$Area=margmeans$Log1P.SunArea.


#rename column for convenience
margmeans$Area=margmeans$Log1P.SunArea.

library(ggplot2)
library(cowplot)
########Asymmetric error bands from manually exponentiated mean +-SE
p0<- ggplot(margmeans, aes(x=Area, y=btmean)) + 
  geom_line(size=2) + theme_classic(base_size=24) + 
  theme(axis.line.y=element_line(size=2, color="black"),
        axis.line.x=element_line(size=2, color="black"))
  #add confidence bands
  p01<-p0+geom_ribbon(aes(ymin=btlo, ymax=bthi), alpha=0.5) #pretty!
p01
ylab<-expression(bold(italic(Crithidia)~count~("cells*"~0.02~mu*L^-1)))
xlabel<-  expression(bold(Log[10]~(Sunflower~area~(m^2))))
p1<-p01 +  theme_cowplot(font_size = 25) +
  xlab(xlabel)+ylab(ylab)
p2<-p1 +  theme(text=element_text(face="bold"), 
        line=element_line(size=2)) +
  theme(axis.line.x=element_line(size=2),
        axis.line.y=element_line(size=2))
p3<- p2 + scale_x_continuous(expand = c(0,0))#removes gap
p3

#now add the raw data
scattered<-p2 + geom_point(data=farmdata,
                aes(x=Log1P.SunArea., y=Count))
scattered

#Export figure:
##for epy:
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Analysis.Individual.Bee.and.Farm/Analyses/Farm Sampling/Farm Sampling Figures")
# ggsave("Farmfig_v4_tmbSE.pdf",height=7, width=7)

################################################################
####update 2017.01
##jg requested plot with tighter y-axis scale
#a log-scaling of y-axis would make the model appear less diffuse
#####################################################################

(lnplot.autoscale<-scattered + scale_y_continuous(trans = "log1p", 
                                                  breaks = c(0, 5, 20, 80, 240)) )


# ggsave("farmplot.logaxis.TMB.pdf", height = 7, width = 7)

#copy to Main figures folder
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Manuscript Drafts/Figures")
# ggsave("farmplot.logaxis.TMB.pdf", height = 7, width = 7)



#######################################################
####################### Appendix #######################                     
#####################################################
##Checks:
#refit with glmer.nb and recreate plot:
  #This is to confirm that lsmeans are anomalous from glmmADMB but not glmmTMB
library(lme4)
m.lme<-glmer.nb(Count ~ Log1P.SunArea. + CollectionDate +(1|Farm), 
                    data=farmdata) 
warnings()
plot(m.lme)

#lsmeans:
lsm.lme<-summary(lsmeans(m.lme, ~Log1P.SunArea., 
                           at=list(Log1P.SunArea.=areas)))
#Manually exponentiate the (mean+-SE) from scale of linear predictor
lsm.lme$btmean<-exp(lsm.lme$lsmean)
lsm.lme$rawhi<-lsm.lme$lsmean + lsm.lme$SE
lsm.lme$rawlo<-lsm.lme$lsmean - lsm.lme$SE
lsm.lme$bthi<-exp(lsm.lme$rawhi)
lsm.lme$btlo<-exp(lsm.lme$rawlo)                        
lsm.lme$Area=lsm.lme$Log1P.SunArea.

plot.lme<-ggplot(lsm.lme, aes(x=Area, y=btmean)) + 
  geom_line(size=2) + theme_classic(base_size=24) + 
  theme(axis.line.y=element_line(size=2, color="black"),
        axis.line.x=element_line(size=2, color="black")) +
  ylab(ylab) + xlab(xlabel)
#add confidence bands
plot.L2<-plot.lme+geom_ribbon(aes(ymin=btlo, ymax=bthi), alpha=0.5) #pretty!
plot.L2

#Compare to glmmTMB fit:
library(cowplot)
plot_grid(p3, plot.L2, align = "h", labels = c("glmmTMB", "glmer.nb"))
#pretty similar
#Export:
# setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Analysis.Individual.Bee.and.Farm/Analyses/Farm Sampling/Farm Sampling Figures")
# ggsave("farmdata.TMB.vs.glmer.nb.pdf", height = 7, width = 10)


citation("glmmTMB")
