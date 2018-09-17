#Script for Jonathan's Consumption Data

library(lme4)
library(MASS)
library(nlme)
library(car)

Pcondata<-read.csv("/Users/rg9403/Desktop/PollenEXP/PollenExp_Winter2014/DATA/CSV/Adj.Net.Pollen.Consump.csv",header=TRUE)
str(Pcondata)
shapiro.test(Pcondata$NetPollenConsumption)
# Shapiro-Wilk normality test
# 
# data:  Pcondata$Adj.Net.Pollen.Consump
# W = 0.9767, p-value = 0.01367 #### Does this matter?

d<- data.frame(Pcondata$NetPollenConsumption,Pcondata$WingSize,Pcondata$CallowMass)
cor(d,use="complete.obs")
#                                   Pcondata.NetPollenConsumption Pcondata.WingSize Pcondata.CallowMass
# Pcondata.NetPollenConsumption                     1.0000000         0.1169731           0.1562572
# Pcondata.WingSize                                 0.1169731         1.0000000           0.9471568
# Pcondata.CallowMass                               0.1562572         0.9471568           1.0000000



########################################################################################################
                    #Linear mixed-model (lme4 package)
########################################################################################################



mod1<-lmer(NetPollenConsumption~Treatment+(1|ColonyID),data=Pcondata)
summary(mod1)
# Linear mixed model fit by REML ['lmerMod']
# Formula: NetPollenConsumption ~ Treatment + (1 | ColonyID)
# Data: Pcondata
# 
# REML criterion at convergence: -1146.3
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -3.5445 -0.5797 -0.0796  0.5757  4.4235 
# 
# Random effects:
#   Groups   Name        Variance  Std.Dev.
# ColonyID (Intercept) 2.089e-06 0.001445
# Residual             7.836e-05 0.008852
# Number of obs: 180, groups:  ColonyID, 8
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)    0.040117   0.001450  27.663
# TreatmentMix   0.054624   0.001869  29.231
# TreatmentRape  0.038128   0.001877  20.311
# TreatmentSun  -0.013734   0.001878  -7.315
# 
# Correlation of Fixed Effects:
#   (Intr) TrtmnM TrtmnR
# TreatmentMx -0.652              
# TreatmentRp -0.653  0.509       
# TreatmentSn -0.651  0.509  0.506

Anova(mod1)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: NetPollenConsumption
# Chisq Df Pr(>Chisq)    
# Treatment 1768.6  3  < 2.2e-16 ***

#Post hoc test for linear mixed models thats works!!!
library(multcomp)
summary(glht(mod1,linfct=mcp(Treatment="Tukey")))

# 
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Tukey Contrasts
# 
# 
# Fit: lmer(formula = NetPollenConsumption ~ Treatment + (1 | ColonyID), 
#           data = Pcondata)
# 
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# Mix - Buck == 0   0.054624   0.001869  29.231   <1e-10 ***
#   Rape - Buck == 0  0.038128   0.001877  20.311   <1e-10 ***
#   Sun - Buck == 0  -0.013734   0.001878  -7.315   <1e-10 ***
#   Rape - Mix == 0  -0.016496   0.001857  -8.884   <1e-10 ***
#   Sun - Mix == 0   -0.068357   0.001857 -36.817   <1e-10 ***
#   Sun - Rape == 0  -0.051862   0.001867 -27.785   <1e-10 ***

AIC(mod1)
# [1] -1134.302

Pcondata1<-na.omit(Pcondata)
mod2<-lmer(NetPollenConsumption~Treatment+WingSize+(1|InocDate)+(1|ColonyID),data=Pcondata1)
summary(mod2)
# Linear mixed model fit by REML ['lmerMod']
# Formula: NetPollenConsumption ~ Treatment + WingSize + (1 | InocDate) +      (1 | ColonyID)
# Data: Pcondata1
# 
# REML criterion at convergence: -1045.6
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -4.1099 -0.5358 -0.0380  0.5733  3.4993 
# 
# Random effects:
#   Groups   Name        Variance  Std.Dev.
# InocDate (Intercept) 9.813e-06 0.003133
# ColonyID (Intercept) 3.003e-07 0.000548
# Residual             6.102e-05 0.007811
# Number of obs: 162, groups:  InocDate, 26; ColonyID, 8
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)    0.006686   0.007161   0.934
# TreatmentMix   0.054578   0.001768  30.874
# TreatmentRape  0.038601   0.001784  21.633
# TreatmentSun  -0.013027   0.001709  -7.622
# WingSize       0.020466   0.004187   4.888
# 
# Correlation of Fixed Effects:
#   (Intr) TrtmnM TrtmnR TrtmnS
# TreatmentMx -0.119                     
# TreatmentRp -0.205  0.473              
# TreatmentSn -0.190  0.484  0.480       
# WingSize    -0.981  0.003  0.093  0.070

Anova(mod2)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: NetPollenConsumption
# Chisq Df Pr(>Chisq)    
# Treatment 1919.261  3  < 2.2e-16 ***
#   WingSize    23.894  1  1.018e-06 ***

AIC(mod2)
#[1] -1029.62

mod3<-lmer(NetPollenConsumption~Treatment+WingSize+Count+(1|InocDate)+(1|ColonyID),data=Pcondata1)
summary(mod3)
# Linear mixed model fit by REML ['lmerMod']
# Formula: Adj.Net.Pollen.Consump ~ Treatment + WingSize + Count + (1 |      InocDate) + (1 | ColonyID)
# Data: Pcondata1
# 
# REML criterion at convergence: -824.1
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -4.2607 -0.6337 -0.0328  0.6885  2.7003 
# 
# Random effects:
#   Groups   Name        Variance  Std.Dev.
# InocDate (Intercept) 1.075e-05 0.003279
# ColonyID (Intercept) 0.000e+00 0.000000
# Residual             5.002e-05 0.007072
# Number of obs: 129, groups:  InocDate, 26; ColonyID, 4
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)   -4.183e-02  8.022e-03  -5.214
# TreatmentMix   1.733e-03  1.984e-03   0.874
# TreatmentRape  1.090e-02  1.953e-03   5.580
# TreatmentSun   8.058e-03  1.963e-03   4.104
# WingSize       2.893e-02  4.603e-03   6.285
# Count          1.920e-05  9.576e-06   2.005
Anova(mod3)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: NetPollenConsumption
# Chisq Df Pr(>Chisq)    
# Treatment 1962.2217  3    < 2e-16 ***
#   WingSize    26.5265  1    2.6e-07 ***
#   Count        4.7812  1    0.02877 * 



AIC(mod2)
# [1] -1029.62
AIC(mod1)
# [1] -1134.302  #### What is a negatove AIC? Good or bad?
AIC(mod3)
#[1] -1010.992


summary(glht(mod3,linfct=mcp(Treatment="Tukey")))
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Tukey Contrasts
# 
# 
# Fit: lmer(formula = Adj.Net.Pollen.Consump ~ Treatment + WingSize + 
#             Count + (1 | InocDate) + (1 | ColonyID), data = Pcondata1)
# 
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# Mix - Buck == 0   0.001733   0.001984   0.874  0.81824    
# Rape - Buck == 0  0.010898   0.001953   5.580  < 0.001 ***
#   Sun - Buck == 0   0.008058   0.001963   4.104  < 0.001 ***
#   Rape - Mix == 0   0.009165   0.001875   4.887  < 0.001 ***
#   Sun - Mix == 0    0.006325   0.001771   3.571  0.00209 ** 
#   Sun - Rape == 0  -0.002841   0.001862  -1.526  0.42151 

########################################################################################################
 #Generalized Linear Mixed-model (glmmADMB package)
########################################################################################################

Pcondata1<-na.omit(Pcondata)
View(Pcondata1) #Looks good.  There are no missing values

library(glmmADMB)

Pcondata1$InocDate<-factor(Pcondata1$InocDate) # all random effects must be a factor

fit1<- glmmadmb(NetPollenConsumption ~ Treatment  + WingSize + (1|ColonyID)  + (1|InocDate), 
                    family = "gaussian",
                    data=Pcondata1) 
summary(fit1)
# Call:
#   glmmadmb(formula = NetPollenConsumption ~ Treatment + WingSize + 
#              (1 | ColonyID) + (1 | InocDate), data = Pcondata1, family = "gaussian")
# 
# AIC: -1084.1 
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    0.00661    0.00710    0.93     0.35    
# TreatmentMix   0.05459    0.00174   31.29  < 2e-16 ***
#   TreatmentRape  0.03859    0.00176   21.90  < 2e-16 ***
#   TreatmentSun  -0.01302    0.00169   -7.72  1.2e-14 ***
#   WingSize       0.02050    0.00415    4.94  7.9e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Number of observations: total=162, ColonyID=8, InocDate=26 
# Random effect variance(s):
#   Group=ColonyID
# Variance    StdDev
# (Intercept) 1.115e-07 0.0003339
# Group=InocDate
# Variance   StdDev
# (Intercept) 9.018e-06 0.003003
# 
# Residual variance: 0.0077138 (std. err.: 0.00046332)

Log-likelihood: 550.073 

summary(glht(fit1,linfct=mcp(Treatment="Tukey")))
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Tukey Contrasts
# 
# 
# Fit: glmmadmb(formula = NetPollenConsumption ~ Treatment + WingSize + 
#                 (1 | ColonyID) + (1 | InocDate), data = Pcondata1, family = "gaussian")
# 
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# Mix - Buck == 0   0.054587   0.001745  31.289  < 1e-10 ***
#   Rape - Buck == 0  0.038590   0.001762  21.896  < 1e-10 ***
#   Sun - Buck == 0  -0.013019   0.001687  -7.715  < 1e-10 ***
#   Rape - Mix == 0  -0.015997   0.002474  -6.467 2.37e-10 ***
#   Sun - Mix == 0   -0.067606   0.002427 -27.850  < 1e-10 ***
#   Sun - Rape == 0  -0.051609   0.002450 -21.067  < 1e-10 ***

fit2<- glmmadmb(NetPollenConsumption ~ Treatment  + WingSize + Count+ (1|ColonyID)  + (1|InocDate), 
                family = "gaussian",
                data=Pcondata1) 
summary(fit2)

# Call:
#   glmmadmb(formula = NetPollenConsumption ~ Treatment + WingSize + 
#              Count + (1 | ColonyID) + (1 | InocDate), data = Pcondata1, 
#            family = "gaussian")
# 
# AIC: -1086.9 
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    2.96e-03   7.19e-03    0.41    0.681    
# TreatmentMix   5.64e-02   1.90e-03   29.67  < 2e-16 ***
#   TreatmentRape  4.00e-02   1.85e-03   21.65  < 2e-16 ***
#   TreatmentSun  -1.10e-02   1.90e-03   -5.78  7.4e-09 ***
#   WingSize       2.15e-02   4.11e-03    5.22  1.8e-07 ***
#   Count          2.04e-05   9.28e-06    2.20    0.028 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Number of observations: total=162, ColonyID=8, InocDate=26 
# Random effect variance(s):
#   Group=ColonyID
# Variance    StdDev
# (Intercept) 1.042e-07 0.0003227
# Group=InocDate
# Variance   StdDev
# (Intercept) 9.753e-06 0.003123
# 
# Residual variance: 0.0075713 (std. err.: 0.00045384)
# 
# Log-likelihood: 552.452 

summary(glht(fit2,linfct=mcp(Treatment="Tukey")))
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Tukey Contrasts
# 
# 
# Fit: glmmadmb(formula = NetPollenConsumption ~ Treatment + WingSize + 
#                 Count + (1 | ColonyID) + (1 | InocDate), data = Pcondata1, 
#               family = "gaussian")
# 
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# Mix - Buck == 0   0.056389   0.001900  29.672  < 1e-08 ***
#   Rape - Buck == 0  0.040003   0.001848  21.650  < 1e-08 ***
#   Sun - Buck == 0  -0.010976   0.001899  -5.781 1.77e-08 ***
#   Rape - Mix == 0  -0.016386   0.002642  -6.202  < 1e-08 ***
#   Sun - Mix == 0   -0.067364   0.002687 -25.074  < 1e-08 ***
#   Sun - Rape == 0  -0.050978   0.002661 -19.158  < 1e-08 ***

library(lsmeans)
My.Lsmeans <- lsmeans(fit2, ~Treatment) 
My.Lsmeans
# Treatment     lsmean          SE df   asymp.LCL  asymp.UCL
# Buck      0.03917529 0.009609762 NA 0.020340503 0.05801008
# Mix       0.09556420 0.009797366 NA 0.076361716 0.11476668
# Rape      0.07917797 0.009808535 NA 0.059953592 0.09840234
# Sun       0.02819979 0.009798252 NA 0.008995565 0.04740401

library(ggplot2)
#make a data frame from the lsmeans object
Lsmeans.df<-summary(My.Lsmeans)
str(Lsmeans.df)
View(Lsmeans.df)
#label for x axis
pollen.axis<-c( "Buckwheat","Mix","Rape", "Sunflower")

plot0<- ggplot(Lsmeans.df, aes(x=Treatment, y=lsmean)) + 
  geom_bar(position=position_dodge(), stat="identity", colour="black" ) +
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE),
                size=1, width=0.2,  # Width of the error bars
                position=position_dodge(0.9))
plot0 #Solid base plot
plot01<-plot0 + theme_bw() #remove colored background
ylabel<- expression(bold(paste("Net Pollen Consumption (g)", sep="") )) #check for encoding
plot02<-plot01 + ylab(ylabel) #y label
plot02
#relabel x axis:
plot05<-plot02 + xlab("Pollen type") + # x axis label
  scale_x_discrete(breaks=levels(Lsmeans.df$Treatment) , labels=pollen.axis)
plot05
#Larger font:
plot06<-plot05 + theme(text=element_text (size=25, face="bold") )
plot06
##remove gridlines (option), remove border
neater<- plot06 +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.line = element_line(size = 1.5, colour = "black"))
neater
#axis scaling
Zeroed<- neater + coord_cartesian(ylim = c(0, 0.1))
Zeroed

#if we want y axis to start at 0:
#this removes space between x axis and bottom of bar (use for non-log scale)
Justified<- neater +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.11))
Justified
