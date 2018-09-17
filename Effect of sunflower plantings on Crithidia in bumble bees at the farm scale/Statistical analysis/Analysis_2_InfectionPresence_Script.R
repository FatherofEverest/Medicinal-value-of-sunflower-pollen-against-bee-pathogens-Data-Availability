
#####################################################################################################
############################### Jan 14th, 2016 ######################################################
#####################################################################################################


data2<-read.csv("/Users/rg9403/Desktop/sunflowersamplepatrickdata/SUNFLOWER_FARM_SAMPLING/Sunflower Sampling Data_CSVs/Analysis_2_InfectionPresence.csv",header=TRUE)
str(data2)
# 'data.frame':	667 obs. of  14 variables:
#   $ CollectionDate   : int  208 208 208 208 208 208 208 208 208 208 ...
# $ Farm             : Factor w/ 22 levels "AST","ATL","BRK",..: 17 17 17 17 17 17 17 17 17 17 ...
# $ BeeID            : int  1 2 3 4 5 6 7 8 9 10 ...
# $ WingSIze         : num  2.4 2.92 2.92 2.54 3.19 ...
# $ Count            : int  0 0 103 0 92 0 219 0 0 0 ...
# $ Log1P.Count.     : num  0 0 2.02 0 1.97 ...
# $ FlowerID2        : Factor w/ 2 levels "NONSUN","SUN": 1 2 2 2 2 2 1 1 1 1 ...
# $ Log1P.SunArea.   : num  1.84 1.84 1.84 1.84 1.84 ...
# $ Sun.Area         : num  67.9 67.9 67.9 67.9 67.9 ...
# $ Log1P.FarmArea.  : num  4.52 4.52 4.52 4.52 4.52 ...
# $ FarmArea         : num  33026 33026 33026 33026 33026 ...
# $ CutFlowerArea    : num  5660 5660 5660 5660 5660 ...
# $ TomatoArea       : num  856 856 856 856 856 ...
# $ InfectionPresence: Factor w/ 2 levels "NO","YES": 1 1 2 1 2 1 2 1 1 1 ...

contrasts(data2$InfectionPresence)
#     YES
# NO    0
# YES   1


######Does collection date corrleate with Sunflower area?###########

d <- data.frame(data2$CollectionDate,data2$Log1P.SunArea.,data2$WingSIze)
cor(d,use="complete.obs")
#                          data2.CollectionDate data2.Log1P.SunArea. data2.WingSIze
# data2.CollectionDate            1.0000000           0.02112470     0.18897770
# data2.Log1P.SunArea.            0.0211247           1.00000000    -0.08091181
# data2.WingSIze                  0.1889777          -0.08091181     1.00000000


#######################################################################################################

library(lme4)
mod1<-glmer(InfectionPresence~Log1P.SunArea.+(1|Farm),family=binomial,data=data2)
summary(mod1)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: InfectionPresence ~ Log1P.SunArea. + (1 | Farm)
# Data: data2
# 
# AIC      BIC   logLik deviance df.resid 
# 854.6    868.1   -424.3    848.6      664 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.1634 -0.7605 -0.4806  0.8476  2.1494 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# Farm   (Intercept) 0.7835   0.8852  
# Number of obs: 667, groups:  Farm, 22
# 
# Fixed effects:
#                 Estimate Std. Error z value Pr(>|z|)
# (Intercept)      0.2448     0.3581   0.684    0.494
# Log1P.SunArea.  -0.2907     0.1855  -1.567    0.117
# 
# Correlation of Fixed Effects:
#              (Intr)
# Log1P.SnAr. -0.813

#### Feb 10th, 2016 Let's see if WingSize matters ####
fit1<-glmer(InfectionPresence~WingSIze+(1|Farm),family=binomial,data=data2)
summary(fit1)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: InfectionPresence ~ WingSIze + (1 | Farm)
# Data: data2
# 
# AIC      BIC   logLik deviance df.resid 
# 836.9    850.4   -415.5    830.9      657 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.6479 -0.7643 -0.4841  0.8316  2.3638 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# Farm   (Intercept) 0.8987   0.948   
# Number of obs: 660, groups:  Farm, 22
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   2.1920     0.9770   2.244   0.0249 *
#   WingSIze     -0.9121     0.3599  -2.535   0.0113 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr)
# WingSIze -0.974


#### WingSize is important so let's keep it in the model. Let's do SunArea and Wing Size ####
fit2<-glmer(InfectionPresence~WingSIze+Log1P.SunArea.+(1|Farm),family=binomial,data=data2)
summary(fit2)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: InfectionPresence ~ WingSIze + Log1P.SunArea. + (1 | Farm)
# Data: data2
# 
# AIC      BIC   logLik deviance df.resid 
# 836.4    854.4   -414.2    828.4      656 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.5050 -0.7596 -0.4834  0.8306  2.3366 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# Farm   (Intercept)    0.7721   0.8787  
# Number of obs: 660, groups:  Farm, 22
# 
# Fixed effects:
#                   Estimate Std. Error z value Pr(>|z|)   
#  (Intercept)      2.7302     1.0251   2.663  0.00774 **
#   WingSIze        -0.9369     0.3595  -2.606  0.00917 **
#   Log1P.SunArea.  -0.3008     0.1848  -1.628  0.10347   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#              (Intr) WngSIz
# WingSIze    -0.938       
# Log1P.SnAr. -0.316  0.036

fit3<-glmer(InfectionPresence~WingSIze+Log1P.SunArea.+Log1P.SunArea.*WingSIze+(1|Farm),family=binomial,data=data2)
summary(fit3)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: InfectionPresence ~ WingSIze + Log1P.SunArea. + Log1P.SunArea. *      WingSIze + (1 | Farm)
# Data: data2
# 
# AIC      BIC   logLik deviance df.resid 
# 833.2    855.7   -411.6    823.2      655 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.4560 -0.7454 -0.4672  0.8462  2.2629 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# Farm   (Intercept) 0.7711   0.8781  
# Number of obs: 660, groups:  Farm, 22
# 
# Fixed effects:
#                            Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                 5.9006     1.7686   3.336 0.000849 ***
#   WingSIze                 -2.1059     0.6406  -3.287 0.001011 ** 
#   Log1P.SunArea.           -2.3692     0.9374  -2.528 0.011487 *  
#   WingSIze:Log1P.SunArea.   0.7676     0.3404   2.255 0.024131 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) WngSIz L1P.SA
# WingSIze    -0.979              
# Log1P.SnAr. -0.831  0.812       
# WSI:L1P.SA.  0.811 -0.824 -0.980












mod2<-glmer(InfectionPresence~Log1P.SunArea.+CollectionDate+WingSIze+(1|Farm),family=binomial,data=data2)
summary(mod2)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: InfectionPresence ~ Log1P.SunArea. + CollectionDate + WingSIze +      (1 | Farm)
# Data: data2
# 
# AIC      BIC   logLik deviance df.resid 
# 831.2    853.7   -410.6    821.2      655 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.4178 -0.7580 -0.5111  0.8393  2.5619 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# Farm   (Intercept) 0.5008   0.7076  
# Number of obs: 660, groups:  Farm, 22
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    11.43134    3.16650   3.610 0.000306 ***
#   Log1P.SunArea. -0.28437    0.15578  -1.825 0.067939 .  
# CollectionDate -0.03871    0.01344  -2.880 0.003978 ** 
#   WingSIze       -0.88946    0.35567  -2.501 0.012391 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) L1P.SA CllctD
# Log1P.SnAr. -0.072              
# CollectinDt -0.949 -0.018       
# WingSIze    -0.237  0.044 -0.065



#####################################################################################################
#####################################################################################################
#####################################################################################################





res=residuals(mod17, type="pearson")

par(mfrow = c(2,2))
plot(fitted(mod17),res)

plot(fitted(mod17), resid(mod17))
