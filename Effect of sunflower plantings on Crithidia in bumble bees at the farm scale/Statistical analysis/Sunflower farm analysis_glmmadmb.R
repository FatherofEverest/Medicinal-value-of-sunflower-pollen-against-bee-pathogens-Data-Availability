farmdata<-read.csv("/Users/rg9403/Desktop/sunflowersamplepatrickdata/SUNFLOWER_FARM_SAMPLING/Sunflower Sampling Data_CSVs/Analysis_1_All.csv",header=TRUE)
str(farmdata)
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
d<- data.frame(farmdata$CollectionDate,farmdata$Log1P.SunArea.,farmdata$WingSize)
cor(d,use="complete.obs")
#                          farmdata.CollectionDate farmdata.Log1P.SunArea. farmdata.WingSize
# farmdata.CollectionDate               1.0000000              0.02159180        0.18897770
# farmdata.Log1P.SunArea.               0.0215918              1.00000000       -0.08091181
# farmdata.WingSize                     0.1889777             -0.08091181        1.00000000


farmdata$CollectionDate<-factor(farmdata$CollectionDate)
farmdata$Count<-as.numeric(farmdata$Count)

Data<-na.omit(farmdata) #in order for glmmadmb's to run

#################################### POISSON #######################################
fitPoiss<- glmmadmb(Count ~ Log1P.SunArea.  + WingSize + (1|Farm)  , 
                    family = "poisson",
                    data=Data)  
summary(fitPoiss)
# 
# Call:
#   glmmadmb(formula = Count ~ Log1P.SunArea. + WingSize + (1 | Farm), 
#            data = Data, family = "poisson")
# 
# AIC: 10406.6 
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)       7.845      0.662    11.8   <2e-16 ***
#   Log1P.SunArea.   -0.470      0.313    -1.5     0.13    
# WingSize         -2.775      0.114   -24.2   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Number of observations: total=660, Farm=22 
# Random effect variance(s):
#   Group=Farm
# Variance StdDev
# (Intercept)    2.476  1.573
# 
# 
# Log-likelihood: -5199.28 


################################ ZERO-INFLATED NEGATIVE BINOMIAL ####################################
zinfl<- glmmadmb(Count ~ Log1P.SunArea.  + WingSize + (1|Farm), 
                 family = "nbinom",zeroInflation = TRUE, data=Data) 

summary(zinfl)
# Call:
#   glmmadmb(formula = Count ~ Log1P.SunArea. + WingSize + (1 | Farm), 
#            data = Data, family = "nbinom", zeroInflation = TRUE)
# 
# AIC: 3796.6 
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      3.8631     0.7691    5.02  5.1e-07 ***
#   Log1P.SunArea.  -0.1836     0.1002   -1.83    0.067 .  
# WingSize         0.0166     0.2844    0.06    0.953    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Number of observations: total=660, Farm=22 
# Random effect variance(s):
#   Group=Farm
# Variance StdDev
# (Intercept)   0.1361 0.3689
# 
# Negative binomial dispersion parameter: 0.57649 (std. err.: 0.072668)
# Zero-inflation: 0.49131  (std. err.:  0.025081 )
# 
# Log-likelihood: -1892.3 



#################################### NEGATIVE BINOMIAL #################################################
nbinom0<-glmmadmb(Count ~ Log1P.SunArea.  + WingSize + (1|Farm), 
                  family = "nbinom",
                  data=Data) 
summary(nbinom0) 

# Call:
#   glmmadmb(formula = Count ~ Log1P.SunArea. + WingSize + (1 | Farm), 
#            data = Data, family = "nbinom")
# 
# AIC: 3843.3 
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)       4.397      1.128    3.90  9.7e-05 ***
#   Log1P.SunArea.   -0.309      0.142   -2.18    0.029 *  
#   WingSize         -0.394      0.410   -0.96    0.337    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Number of observations: total=660, Farm=22 
# Random effect variance(s):
#   Group=Farm
# Variance StdDev
# (Intercept)   0.2814 0.5305
# 
# Negative binomial dispersion parameter: 0.13786 (std. err.: 0.0096981)
# 
# Log-likelihood: -1916.67 

m1<-update(nbinom0, ~. - WingSize) 
summary(m1)
# 
# Call:
#   glmmadmb(formula = Count ~ Log1P.SunArea. + (1 | Farm), data = Data, 
#            family = "nbinom")
# 
# AIC: 3842.3 
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)       3.351      0.279   12.02   <2e-16 ***
#   Log1P.SunArea.   -0.308      0.146   -2.11    0.035 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Number of observations: total=660, Farm=22 
# Random effect variance(s):
#   Group=Farm
# Variance StdDev
# (Intercept)   0.3099 0.5567
# 
# Negative binomial dispersion parameter: 0.13791 (std. err.: 0.0096856)
# 
# Log-likelihood: -1917.13 
anova(nbinom0, m1) #if p>0.05, you are good to drop the parameter (in this case WingSize)
####and keep the simpler model
# Analysis of Deviance Table
# 
# Model 1: Count ~ Log1P.SunArea.
# Model 2: Count ~ Log1P.SunArea. + WingSize
# NoPar  LogLik Df Deviance Pr(>Chi)
# 1     4 -1917.1                     
# 2     5 -1916.7  1     0.92   0.3375


##### Since Wing Size doesn't have an effect, let's use the original data because there were a bunch of rows that 
 # we needed to omit due to Wing Size NAs
nbinom2<-glmmadmb(Count ~ Log1P.SunArea. + (1|Farm), 
                  family = "nbinom",
                  data=farmdata) 
summary(nbinom2) 
# Call:
#   glmmadmb(formula = Count ~ Log1P.SunArea. + (1 | Farm), data = farmdata, 
#            family = "nbinom")
# 
# AIC: 3908 
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)       3.355      0.271   12.39   <2e-16 ***
#   Log1P.SunArea.   -0.294      0.142   -2.07    0.038 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Number of observations: total=667, Farm=22 
# Random effect variance(s):
#   Group=Farm
# Variance StdDev
# (Intercept)   0.2808 0.5299
# 
# Negative binomial dispersion parameter: 0.13856 (std. err.: 0.0096654)
# 
# Log-likelihood: -1949.98 

AIC(nbinom0) 
#[1] 3843.34
AIC(m1)
#[1] 3842.26
AIC(nbinom2)
#[1] 3907.96 # This is the model that doesn't omit any NA's because there are no longer any other explanatories

plot(Data$Log1P.SunArea.,Data$Count)



