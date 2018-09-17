#######SCRIPT TO ANALYZE CRITHIDIA FROM WILD BEES########
###Collected from Amherst area farms summer 2015
####PROPOSED EXPLANATORY VARIABLE: sunflower area#

#Updated 2017.05.26 -- lsmeans behaves strangely; feed glmmTMB model instead

rm(list=ls())

library(glmmADMB)
library(glmmTMB)
library(lsmeans)
library(ggplot2)

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
#### Check correlations between explanatory variables, especially collectiondate.####
d<- data.frame(farmdata$CollectionDate,farmdata$Log1P.SunArea.,farmdata$WingSize)
cor(d,use="complete.obs")
#                          farmdata.CollectionDate farmdata.Log1P.SunArea. farmdata.WingSize
# farmdata.CollectionDate               1.0000000              0.02159180        0.18897770
# farmdata.Log1P.SunArea.               0.0215918              1.00000000       -0.08091181
# farmdata.WingSize                     0.1889777             -0.08091181        1.00000000
#none of the variables strongly correlated

#farmdata$CollectionDate<-factor(farmdata$CollectionDate)
#EPY note-- collection date I will keep as numeric covariate....
#to account for increasing/decreasing infection across the sampling period
farmdata$Count<-as.numeric(farmdata$Count)

Data<-na.omit(farmdata) #in order for glmmadmb's to run

Data1<-Data

#simple boxplot
library(car)
Boxplot(Count~Log1P.SunArea., data = Data1)


###Scout some distributions, use Julia Pilowky script
#http://www.juliapilowsky.com/mixedmodels/
#normal
library(MASS)
library(car)
qqp(Data1$Count, "norm") #overdispersed
# lnorm means lognormal
qqp(Data1$Count, "lnorm") 
#for this we would have to transform zeroes to ones

#Poisson:
poisson <- fitdistr(Data1$Count, "Poisson")
qqp(Data1$Count, "pois", poisson$estimate)
#not good, way over-dispersed

#negative binomial
nbinom <- fitdistr(Data1$Count, "Negative Binomial")
qqp(Data1$Count, "nbinom", size = nbinom$estimate[[1]],
    mu = nbinom$estimate[[2]])
#that's pretty terrible, under-dispersed
nbinom <- fitdistr(1+Data1$Count, "Negative Binomial")
qqp(1+Data1$Count, "nbinom", size = nbinom$estimate[[1]],
    mu = nbinom$estimate[[2]]) #worked better

plot(Data$Count~Data$SunArea)
byarea<-lm(Data$Count~Data$SunArea)
abline(byarea)
summary(byarea) 
Anova(byarea) #interesting trend there
#but it is not accounting for the non-independence
#of bees from same farm

plot(Data$Count~Data$Log1P.SunArea.)
bylog<-lm(Data$Count~Data$Log1P.SunArea.) #steep!
abline(bylog)
summary(bylog) #interesting, quite significant trend
#but does not account for non-independence 
#of bees from same farm

#here let's also scout the effects of collection date
plot(Data$CollectionDate,Data$Count)
#decreasing crithidia over time
bydate<-lm(Data$Count~Data$CollectionDate)
abline(bydate)
Anova(bydate)


################################ ZERO-INFLATED NEGATIVE BINOMIAL ####################################
library(glmmADMB)
Log0<- glmmadmb(Count ~ Log1P.SunArea.  + WingSize + 
                   CollectionDate +
                   (1|Farm), 
                 family = "nbinom",zeroInflation = FALSE, data=Data) 
summary(Log0)
#Compare to zero-inflated:
Log.zi<-glmmadmb(Count ~ Log1P.SunArea.  + WingSize + 
                   CollectionDate +
                   (1|Farm), 
                 family = "nbinom",zeroInflation = TRUE, data=Data) 

library(bbmle)
AICtab(Log0, Log.zi) #this is quite a big improvement with zero inflation
#Tough call whether there are 2 processes at work or not
#In some sense, every bee ought to have been exposed to some Crithidia 
  #by the time it starts to forage
#However, one could envision that some bees have not had sufficient time or exposure
   #to develop detectable infection
plot(Log0$residuals)
plot(Log0$fitted, Log0$residuals)
qqp(Log0$residuals) #oof

plot(Log.zi$fitted, Log.zi$residuals)
plot(Log.zi$residuals)
qqp(Log.zi$residuals) #this model also veers off the line at high values

summary(Log0)
summary(Log.zi) #Zero-inflation: 0.4862  (std. err.:  0.026077 )

#Let's drop wing size p=0.6
L1<-update(Log0, ~. -WingSize)
anova(Log0, L1)
AICtab(Log0, L1) #keep the simpler model
str(Data$WingSize)

#now we can reinstate the bees
#for whom we did not have wing measures
L1plus<-update(L1, data=farmdata)
#looks good-- all predictors are helping us there
summary(L1plus)
Anova(L1plus)
# Analysis of Deviance Table (Type II tests)
# Response: Count
# Df   Chisq Pr(>Chisq)    
# Log1P.SunArea.  1  7.3162   0.006833 ** 
#   CollectionDate  1 16.4258  5.059e-05 ***
#so this is the model of choice:
#L1plus<-glmmadmb(Count ~ Log1P.SunArea. + CollectionDate +(1|Farm), 
        #family = "nbinom",zeroInflation = FALSE, data=farmdata) 

#############PLOTTING##############

#REFIT the model with glmmTMB: seems that this is only way to get logical lsmeans standard errors
?glmmTMB
#First check ZI model:
L1plus.tmb.zi<-glmmTMB(Count ~ Log1P.SunArea. + CollectionDate +(1|Farm), 
                    family = "nbinom2",ziformula = ~1, data=farmdata) 
drop1(L1plus.tmb.zi, test = "Chisq")
# Df    AIC    LRT Pr(>Chi)   
# <none>            3850.3                   
# Log1P.SunArea.  1 3852.3 3.9819 0.045991 * 
#   CollectionDate  1 3857.1 8.8672 0.002903 **

#Also the binomial model:
L1plus.tmb.bi<-glmmTMB(Count>0 ~ Log1P.SunArea. + CollectionDate +(1|Farm), 
                       family = "binomial",ziformula = ~0, data=farmdata) 
drop1(L1plus.tmb.bi, test = "Chisq")
# Count > 0 ~ Log1P.SunArea. + CollectionDate + (1 | Farm)
# Df    AIC    LRT Pr(>Chi)   
# <none>            848.94                   
# Log1P.SunArea.  1 849.87 2.9238 0.087284 . 
# CollectionDate  1 854.58 7.6319 0.005735 **


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
