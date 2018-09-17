data<-read.csv("/Users/rg9403/Desktop/sunflowersamplepatrickdata/SUNFLOWER_FARM_SAMPLING/Sunflower Sampling Data_CSVs/Analysis_1_All.csv",header=TRUE)
str(data)

# 'data.frame':	667 obs. of  13 variables:
#   $ CollectionDate : int  208 208 208 208 208 208 208 208 208 208 ...
# $ Farm           : Factor w/ 22 levels "AST","ATL","BRK",..: 17 17 17 17 17 17 17 17 17 17 ...
# $ BeeID          : int  1 2 3 4 5 6 7 8 9 10 ...
# $ WingSize       : num  2.4 2.92 2.92 2.54 3.19 ...
# $ FlowerID       : Factor w/ 18 levels "BAS","BRG","BSUN",..: 2 3 3 3 3 3 7 9 9 9 ...
# $ Count          : int  0 0 103 0 92 0 219 0 0 0 ...
# $ Log1P.Count.   : num  0 0 2.02 0 1.97 ...
# $ FlowerID2      : Factor w/ 2 levels "NONSUN","SUN": 1 2 2 2 2 2 1 1 1 1 ...
# $ FlowerID3      : Factor w/ 3 levels ".","HELI","NONHELI": 3 2 2 2 2 2 3 3 3 3 ...
# $ FarmSize       : Factor w/ 22 levels "112,137.41","131,041.50",..: 10 10 10 10 10 10 10 10 10 10 ...
# $ Log1P.FarmSize.: num  4.52 4.52 4.52 4.52 4.52 4.52 4.52 4.52 4.52 4.52 ...
# $ SunArea        : num  67.9 67.9 67.9 67.9 67.9 ...
# $ Log1P.SunArea. : num  1.84 1.84 1.84 1.84 1.84 ...

#### Check correlations between explanatory variables, especially collectiondate.####
d<- data.frame(data$CollectionDate,data$Log1P.SunArea.,data$WingSize)
cor(d,use="complete.obs")
#                         data.CollectionDate data.Log1P.SunArea. data.WingSize
# data.CollectionDate           1.0000000          0.02159180    0.18897770
# data.Log1P.SunArea.           0.0215918          1.00000000   -0.08091181
# data.WingSize                 0.1889777         -0.08091181    1.00000000
cov(d,use="complete.obs")
#                       data.CollectionDate data.Log1P.SunArea. data.WingSize
# data.CollectionDate         169.2237458          0.31991981    0.67160330
# data.Log1P.SunArea.           0.3199198          1.29730898   -0.02517707
# data.WingSize                 0.6716033         -0.02517707    0.07463512

###For a poisson model, the dependant varibale must be an integer

####Normal Probability Plot of Residuals####
Count.lm = lm(Count ~ Sun.Area+TomatoArea+CutFlowerArea, data=data) 
Count.stdres = rstandard(Count.lm)
qqnorm(Count.stdres,ylab="Standardized Residuals",xlab="Normal Scores",main="Crithidia Count") 
qqline(Count.stdres)
####Short Tails - An S shaped-curve indicates shorter than normal tails, i.e. less variance 
####than expected.


library(lme4)
Count_mod1<-glmer(Count ~ Log1P.SunArea.+(1|BeeID)+(1|Farm),family=poisson(link=log),na.action=na.exclude,data=data)
summary(Count_mod1)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: poisson  ( log )
# Formula: Count ~ Log1P.SunArea. + (1 | BeeID) + (1 | Farm)
# Data: data
# 
# AIC      BIC   logLik deviance df.resid 
# 4043.0   4061.1  -2017.5   4035.0      663 
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -0.59101 -0.38570 -0.27404  0.03643  0.18407 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# BeeID  (Intercept) 11.227   3.351   
# Farm   (Intercept)  1.426   1.194   
# Number of obs: 667, groups:  BeeID, 667; Farm, 22
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)      0.3346     0.5212   0.642   0.5210  
# Log1P.SunArea.  -0.6571     0.2683  -2.449   0.0143 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr)
# Log1P.SnAr. -0.775

Count_mod2<-glmer(Count ~ WingSize+(1|BeeID)+(1|Farm),family=poisson(link=log),data=data)
summary(Count_mod2)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: poisson  ( log )
# Formula: Count ~ WingSize + (1 | BeeID) + (1 | Farm)
# Data: data
# 
# AIC      BIC   logLik deviance df.resid 
# 3968.2   3986.1  -1980.1   3960.2      656 
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -0.60899 -0.38203 -0.27813  0.03519  0.17281 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# BeeID  (Intercept) 11.001   3.317   
# Farm   (Intercept)  1.767   1.329   
# Number of obs: 660, groups:  BeeID, 660; Farm, 22
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)   4.2193     1.7364   2.430  0.01510 * 
#   WingSize     -1.8708     0.6511  -2.873  0.00406 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr)
# WingSize -0.979

Count_mod3<-glmer(Count ~ Log1P.SunArea.+WingSize+(1|BeeID)+(1|Farm),family=poisson(link=log),data=data)
summary(Count_mod3)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: poisson  ( log )
# Formula: Count ~ Log1P.SunArea. + WingSize + (1 | BeeID) + (1 | Farm)
# Data: data
# 
# AIC      BIC   logLik deviance df.resid 
# 3963.9   3986.4  -1977.0   3953.9      655 
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -0.60201 -0.37151 -0.27023  0.03593  0.18605 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# BeeID  (Intercept) 11.264   3.356   
# Farm   (Intercept)  1.186   1.089   
# Number of obs: 660, groups:  BeeID, 660; Farm, 22
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)      5.5921     1.7968   3.112  0.00186 **
#   Log1P.SunArea.  -0.6749     0.2520  -2.678  0.00740 **
#   WingSize        -1.9960     0.6583  -3.032  0.00243 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) L1P.SA
# Log1P.SnAr. -0.243       
# WingSize    -0.962  0.035

Count_mod4<-glmer(Count ~ Log1P.SunArea.+WingSize+Log1P.SunArea.*WingSize+(1|BeeID)+(1|Farm),family=poisson(link=log),data=data)
summary(Count_mod4)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: poisson  ( log )
# Formula: Count ~ Log1P.SunArea. + WingSize + Log1P.SunArea. * WingSize +      (1 | BeeID) + (1 | Farm)
# Data: data
# 
# AIC      BIC   logLik deviance df.resid 
# 3961.5   3988.4  -1974.7   3949.5      654 
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -0.62051 -0.35884 -0.26211  0.03479  0.18405 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# BeeID  (Intercept) 11.188   3.345   
# Farm   (Intercept)  1.179   1.086   
# Number of obs: 660, groups:  BeeID, 660; Farm, 22
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)              10.3555     2.9029   3.567 0.000361 ***
#   Log1P.SunArea.           -4.0146     1.6116  -2.491 0.012735 *  
#   WingSize                 -3.7693     1.0778  -3.497 0.000470 ***
#   Log1P.SunArea.:WingSize   1.2510     0.5951   2.102 0.035558 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#               (Intr) Lg1P.SA. WingSz
# Log1P.SnAr. -0.799                
# WingSize    -0.986  0.785         
# Lg1P.SA.:WS  0.786 -0.988   -0.792



NB_mod1<- glmer.nb(Count ~ Log1P.SunArea.+(1|Farm),data=data)
summary(NB_mod1)

NB_mod2<- glmer.nb(Count ~ WingSize+(1|Farm),data=data)
summary(NB_mod2)

NB_mod3<- glmer.nb(Count ~ Log1P.SunArea.+WingSize+(1|Farm),data=data)
summary(NB_mod3)

# 
# > NB_mod1<- glmer.nb(Count ~ Log1P.SunArea.+(1|Farm),data=data)
# > summary(NB_mod1)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: Negative Binomial(0.1382)  ( log )
# Formula: Count ~ Log1P.SunArea. + (1 | Farm)
# Data: ..2
# 
# AIC      BIC   logLik deviance df.resid 
# 3908.6   3926.6  -1950.3   3900.6      663 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -0.4899 -0.4886 -0.4864  0.1323  7.8738 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Farm     (Intercept) 0.1442   0.3798  
# Residual             0.5745   0.7580  
# Number of obs: 667, groups:  Farm, 22
# 
# Fixed effects:
#   Estimate Std. Error t value Pr(>|z|)    
# (Intercept)      3.2947     0.2643  12.467   <2e-16 ***
#   Log1P.SunArea.  -0.2899     0.1362  -2.128   0.0333 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr)
# Log1P.SnAr. -0.768
# > NB_mod2<- glmer.nb(Count ~ WingSize+(1|Farm),data=data)
# > summary(NB_mod2)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: Negative Binomial(0.1382)  ( log )
# Formula: Count ~ WingSize + (1 | Farm)
# Data: ..2
# 
# AIC      BIC   logLik deviance df.resid 
# 3846.5   3864.5  -1919.3   3838.5      656 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -0.4982 -0.4966 -0.4928  0.1069  6.4500 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Farm     (Intercept) 0.2379   0.4877  
# Residual             0.5554   0.7453  
# Number of obs: 660, groups:  Farm, 22
# 
# Fixed effects:
#   Estimate Std. Error t value Pr(>|z|)    
# (Intercept)   3.7036     1.1078   3.343 0.000828 ***
#   WingSize     -0.3491     0.4117  -0.848 0.396456    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr)
# WingSize -0.985
# > NB_mod3<- glmer.nb(Count ~ Log1P.SunArea.+WingSize+(1|Farm),data=data)
# > summary(NB_mod3)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: Negative Binomial(0.1375)  ( log )
# Formula: Count ~ Log1P.SunArea. + WingSize + (1 | Farm)
# Data: ..2
# 
# AIC      BIC   logLik deviance df.resid 
# 3844.0   3866.4  -1917.0   3834.0      655 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -0.4832 -0.4818 -0.4795  0.1176  8.2295 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Farm     (Intercept) 0.1480   0.3847  
# Residual             0.5876   0.7666  
# Number of obs: 660, groups:  Farm, 22
# 
# Fixed effects:
#   Estimate Std. Error t value Pr(>|z|)    
# (Intercept)      4.3610     1.1227   3.884 0.000103 ***
#   Log1P.SunArea.  -0.3053     0.1366  -2.234 0.025458 *  
#   WingSize        -0.4034     0.4075  -0.990 0.322200    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) L1P.SA
# Log1P.SnAr. -0.207       
# WingSize    -0.972  0.026
