####analysis of Exp 9_ USA_WI_SUNFOWER POLLEN EXPERIMENT ####### 
###individual bees.... JONATHAN GIACOMINI.... 4 POLLEN DIETS EFFECTS ON CRITHIDIA COUNTS#
#1 - SUN - WI Sunflower Pollen
#2 - CSUN - China Sunflower Pollen
#3 - MIX- WI Mix POllen
#4 - KMIX - Koppert Mix Pollen 

####Final model -> nbinom2 (glmmadmb) or M2 (lme4) # not sure which is better


data<-read.csv("/Users/rg9403/Google Drive/Sunflower Pollen Experiments/Irwin Lab_Sunflower Pollen Experiments/Exp 9_ U.S.A._WI_Sunflower Pollen Exp/U.S.A._WI_Sunflower Pollen Exp.csv",header=TRUE)
str(data)
# 'data.frame':	120 obs. of  4 variables:
# $ Bee.ID         : int  1 2 3 4 5 6 7 8 9 10 ...
# $ Colony.ID      : Factor w/ 3 levels "A","B","C": 1 2 3 1 2 3 1 2 3 1 ...
# $ Treatment      : Factor w/ 4 levels "CSUN","KMIX",..: 1 4 3 2 1 4 3 2 1 4 ...
# $ Crithidia.Count: int  39 NA 66 127 NA 0 164 NA 6 NA ...


####Renaming to stay consisitent with Evan's script
Data<-data #to rename 
Data$wingmm<-Data$BeeSize # renames for compatibility
Data$Pollen<-Data$Treatment #ditto
Data$totcrith<-Data$CrithidiaCount #ditto


library(glmmADMB)
Data1<-na.omit(Data) #The following models don't work unless I omit the NAs
#Compare Poisson with zero inflation, regular negative binomial, and negative binomial with zero inflation
#Example:
#fit_zipoiss <- glmmadmb(NCalls ~ (FoodTreatment + ArrivalTime) * SexParent +
#offset(log(BroodSize)) + (1 | Nest), data = Owls, zeroInflation = TRUE,
#family = "poisson") 
#first try poisson model:
fitPoiss<- glmmadmb(totcrith ~ Pollen   + wingmm + (1|ColonyID), 
                    family = "poisson",
                    data=Data1)  
summary(fitPoiss) # Bee Size (wingmm) not significant p = 0.1178

#lets re-run without wingmm and use Data to bring back the one excluded NA from the missing wing measurement
#need to make new Data2 to exclude wingmm (BeeSize) so that we can remove NA's from dead bees that are missing Crithidia Counts
str(data)
data$BeeSize<-NULL
str(data) # ok...BeeSize removed
# lets omit NA's from CrithidiaCount
data2<-na.omit(data)
str(data2)

####Renaming to stay consisitent with Evan's script
Data2<-data2 #to rename 
Data2$wingmm<-Data2$BeeSize # renames for compatibility
Data2$Pollen<-Data2$Treatment #ditto
Data2$totcrith<-Data2$CrithidiaCount #ditto
str(Data2)

fitPoiss2<- glmmadmb(totcrith ~ Pollen   + (1|ColonyID), 
                    family = "poisson",
                    data=Data2) 
summary(fitPoiss2)
# Call:
# glmmadmb(formula = totcrith ~ Pollen + (1 | Colony.ID), data = Data2, 
#          family = "poisson")
# 
# AIC: 2229.6 
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    1.348      0.861    1.57    0.117    
# PollenKMIX     2.044      0.110   18.50   <2e-16 ***
#   PollenMIX      1.441      0.108   13.37   <2e-16 ***
#   PollenSUN      0.326      0.154    2.12    0.034 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Number of observations: total=110, Colony.ID=3 
# Random effect variance(s):
#   Group=Colony.ID
# Variance StdDev
# (Intercept)    2.156  1.468
# 
# 
# Log-likelihood: -1109.78 

library(multcomp)
summary(glht(fitPoiss2,linfct=mcp(Pollen="Tukey")))   
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# KMIX - CSUN == 0  2.04365    0.11048  18.498   <0.001 ***
#   MIX - CSUN == 0   1.44121    0.10780  13.369   <0.001 ***
#   SUN - CSUN == 0   0.32648    0.15373   2.124    0.121    
# MIX - KMIX == 0  -0.60244    0.03656 -16.480   <0.001 ***
#   SUN - KMIX == 0  -1.71717    0.18031  -9.523   <0.001 ***
#   SUN - MIX == 0   -1.11474    0.16685  -6.681   <0.001 ***



#######################################  NEG BIN ########################################################
#####now negative binomial
#Start with full nesting
#Start with Data1 that excludes missing wingmm and missing Crithidia Counts
nbinom0<-glmmadmb(totcrith ~ Pollen   + wingmm + (1|ColonyID)  , 
                  family = "nbinom",
                  data=Data1)  
summary(nbinom0) #wingmm no effect p = 0.94499
# Call:
#   glmmadmb(formula = totcrith ~ Pollen + wingmm + (1 | ColonyID), 
#            data = Data1, family = "nbinom")
# 
# AIC: 988.8 
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   2.4414     1.8263    1.34  0.18129    
# PollenKMIX    1.9557     0.3937    4.97  6.8e-07 ***
#   PollenMIX     1.4676     0.3813    3.85  0.00012 ***
#   PollenSUN     0.2343     0.4096    0.57  0.56739    
# wingmm        0.0399     0.5777    0.07  0.94499    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Number of observations: total=109, ColonyID=3 
# Random effect variance(s):
#   Group=ColonyID
# Variance StdDev
# (Intercept)   0.3262 0.5712
# 
# Negative binomial dispersion parameter: 0.5121 (std. err.: 0.073507)
# 
# Log-likelihood: -487.393 

nbinom2<-glmmadmb(totcrith ~ Pollen + (1|ColonyID)  , 
                  family = "nbinom",
                  data=Data2)  
summary(nbinom2)

# Call:
#   glmmadmb(formula = totcrith ~ Pollen + (1 | Colony.ID), data = Data2, 
#            family = "nbinom")
# 
# AIC: 999.1 
# 
# Coefficients:
#                Estimate  Std. Error z value Pr(>|z|)    
#  (Intercept)    2.564      0.429    5.98  2.2e-09 ***
#   PollenKMIX     1.938      0.378    5.13  2.9e-07 ***
#   PollenMIX      1.466      0.378    3.88   0.0001 ***
#   PollenSUN      0.233      0.406    0.57   0.5660    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Number of observations: total=110, Colony.ID=3 
# Random effect variance(s):
#   Group=Colony.ID
# Variance StdDev
# (Intercept)   0.3254 0.5704
# 
# Negative binomial dispersion parameter: 0.51832 (std. err.: 0.074036)
# 
# Log-likelihood: -493.541 

library(multcomp)
#post hoc test 
summary(glht(nbinom2,linfct=mcp(Pollen="Tukey")))
# 
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Tukey Contrasts
# 
# 
# Fit: glmmadmb(formula = totcrith ~ Pollen + (1 | Colony.ID), data = Data2, 
#               family = "nbinom")
# 
# Linear Hypotheses:
#                   Estimate Std. Error z value Pr(>|z|)    
# KMIX - CSUN == 0   1.9378     0.3778   5.129   <0.001 ***
# MIX - CSUN == 0    1.4660     0.3776   3.882   <0.001 ***
# SUN - CSUN == 0    0.2329     0.4057   0.574   0.9331    
# MIX - KMIX == 0   -0.4718     0.5291  -0.892   0.7918    
# SUN - KMIX == 0   -1.7049     0.5641  -3.022   0.0116 *  
# SUN - MIX == 0    -1.2331     0.5609  -2.198   0.1126    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)



library(car)
Anova(nbinom2) #this will give chisquare tests, similar to an anova table, called "Deviance Table"
# Analysis of Deviance Table (Type II tests)
# 
# Response: totcrith
# Df  Chisq Pr(>Chisq)    
# Pollen      3 41.368  5.464e-09 ***
#   Residuals 104                      



###try with a different package for comparison
library(lme4)
M0<- glmer.nb(totcrith ~ Pollen  + wingmm + (1|ColonyID)  , data=Data1)

summary(M0) #Same thing...wingmm not significnat p = 0.945690

Anova(M0)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: totcrith
# Chisq Df Pr(>Chisq)    
# Pollen 34.2150  3  1.785e-07 ***
#   wingmm  0.0046  1     0.9457

# Let's run it again without wingmm
M2<- glmer.nb(totcrith ~ Pollen  + (1|ColonyID)  , data=Data2)

summary(M2)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: Negative Binomial(0.5186)  ( log )
# Formula: totcrith ~ Pollen + (1 | ColonyID)
# Data: ..2
# 
# AIC      BIC   logLik deviance df.resid 
# 999.1   1015.3   -493.6    987.1      104 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -0.8065 -0.7694 -0.2877  0.3558  4.1655 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# ColonyID (Intercept) 0.2592   0.5091  
# Residual             0.7879   0.8877  
# Number of obs: 110, groups:  ColonyID, 3
# 
# Fixed effects:
#   Estimate Std. Error t value Pr(>|z|)    
# (Intercept)   2.5412     0.4295   5.916 3.29e-09 ***
#   PollenKMIX    1.9368     0.3775   5.131 2.88e-07 ***
#   PollenMIX     1.4650     0.3773   3.883 0.000103 ***
#   PollenSUN     0.2321     0.4054   0.573 0.566932    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) PlKMIX PllMIX
# PollenKMIX -0.463              
# PollenMIX  -0.459  0.520       
# PollenSUN  -0.426  0.480  0.495
Anova(M2)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: totcrith
# Chisq Df Pr(>Chisq)    
# Pollen 36.033  3  7.369e-08 ***
summary(glht(M2,linfct=mcp(Pollen="Tukey")))
# Linear Hypotheses:
#                   Estimate  Std.Error z value Pr(>|z|)    
# KMIX - CSUN == 0   1.9368     0.3775   5.131  < 0.001 ***
# MIX - CSUN == 0    1.4650     0.3773   3.883  < 0.001 ***
# SUN - CSUN == 0    0.2321     0.4054   0.573  0.94019    
# MIX - KMIX == 0   -0.4718     0.3698  -1.276  0.57799    
# SUN - KMIX == 0   -1.7047     0.4000  -4.262  < 0.001 ***
# SUN - MIX == 0    -1.2329     0.3940  -3.129  0.00933 ** 

###################################### PLOTS ####################################################
#Use lsmeans to get mean counts
library(lsmeans)
My.Lsmeans <- lsmeans(nbinom2, ~Pollen) 
My.Lsmeans
altLsmeans<- lsmeans(M2, ~Pollen, data=Data2) #substitute the name of your favorite model
#choose the simplest model you can without excluding significant terms
#obviously error messages may narrow down your choices:)
altLsmeans
### Funny how the means are slightly different between nbinom2 (glmmadmb) and M2 (lme4)

################################ My. Lsmeans (nbinom2 model) Plot ###############################
#BELOW adapted FROM CONROY EXPERIMENT
#######PLOTTING : PARASITE LOAD#############
#plotting marginal means:
library(ggplot2)
#make a data frame from the lsmeans object
Lsmeans.df<-summary(My.Lsmeans)
str(Lsmeans.df)
View(Lsmeans.df)
#label for x axis
pollen.axis<-c( "China Sun","Koppert Mix","US Mix", "US Sun")

plot0<- ggplot(Lsmeans.df, aes(x=Pollen, y=lsmean)) + 
  geom_bar(position=position_dodge(), stat="identity", colour="black" ) +
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE),
                size=1, width=0.2,  # Width of the error bars
                position=position_dodge(0.9))
plot0 #Solid base plot
plot01<-plot0 + theme_bw() #remove colored background
ylabel<- expression(bold(paste("Parasite load (ln (cells * 0.02 ", ?L^-1, "))", sep="") )) #check for encoding
plot02<-plot01 + ylab(ylabel) #y label
plot02
#relabel x axis:
plot05<-plot02 + xlab("Pollen type") + # x axis label
  scale_x_discrete(breaks=levels(Lsmeans.df$Pollen) , labels=pollen.axis)
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
Zeroed<- neater + coord_cartesian(ylim = c(0, 5.5))
Zeroed



#if we want y axis to start at 0:
#this removes space between x axis and bottom of bar (use for non-log scale)
Justified<- neater +
  scale_y_continuous(expand = c(0,0), limits = c(0,5.5))
Justified

fig1<-Justified


#citation info

#install.packages("bibtex")
library(bibtex)
write.bib(entry="ggplot2", file = "ggplot2.bib", append = FALSE, verbose = TRUE)
write.bib(entry="coxme", file = "coxme.bib", append = FALSE, verbose = TRUE)
write.bib(entry="survival", file = "survival.bib", append = FALSE, verbose = TRUE)
write.bib(entry="car", file = "car.bib", append = FALSE, verbose = TRUE)
write.bib(entry="lme4", file = "lme4.bib", append = FALSE, verbose = TRUE)
write.bib(entry="MASS", file = "MASS.bib", append = FALSE, verbose = TRUE)
write.bib(entry="lsmeans", file = "lsmeans.bib", append = FALSE, verbose = TRUE)
write.bib(entry="glmmADMB", file = "glmmADMB.bib", append = FALSE, verbose = TRUE)


#################### DATA SUMMARY ############################
##############################################################

summary(Data2)
# Bee.ID       Colony.ID Treatment Crithidia.Count   Pollen      totcrith     
# Min.   :  1.00   A:37      CSUN:27   Min.   :  0.00   CSUN:27   Min.   :  0.00  
# 1st Qu.: 32.25   B:32      KMIX:29   1st Qu.:  2.00   KMIX:29   1st Qu.:  2.00  
# Median : 60.50   C:41      MIX :29   Median : 22.00   MIX :29   Median : 22.00  
# Mean   : 61.49             SUN :25   Mean   : 51.83   SUN :25   Mean   : 51.83  
# 3rd Qu.: 90.75                       3rd Qu.: 85.25             3rd Qu.: 85.25  
# Max.   :120.00                       Max.   :350.00             Max.   :350.00  

st.err <- function(CrithidiaCount) {sd(CrithidiaCount)/sqrt(length(CrithidiaCount))}
MEANS<- aggregate(CrithidiaCount ~ Treatment, data = Data2, mean)
SE <- aggregate(CrithidiaCount ~ Treatment, data = Data2, st.err)
names(MEANS)[names(MEANS)=="CrithidiaCount"] <- "Mean Crithidia Load"
names(SE)[names(SE)=="CrithidiaCount"] <- "SE"
DESCRIPT <- cbind(MEANS,SE)
DESCRIPT # Columns need names
#     Treatment       Mean Crithidia         SE
#       CSUN            14.81481          3.714159
#       KMIX            95.06897          13.888093
#        MIX            68.89655          14.607453
#        SUN            21.84000           5.759537



########################## altlsmeans (M2 model) Plot ######################################
Lsmeans.df1<-summary(altLsmeans)
str(Lsmeans.df1)
View(Lsmeans.df1)
#label for x axis
pollen.axis<-c( "China Sun","Koppert Mix","US Mix", "US Sun")

plotA<- ggplot(Lsmeans.df1, aes(x=Pollen, y=lsmean)) + 
  geom_bar(position=position_dodge(), stat="identity", colour="black" ) +
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE),
                size=1, width=0.2,  # Width of the error bars
                position=position_dodge(0.9))
plotA #Solid base plot
plotB<-plotA + theme_bw() #remove colored background
ylabel<- expression(bold(paste("Parasite load (ln (cells * 0.02 ", ?L^-1, "))", sep="") )) #check for encoding
plotC<-plotB + ylab(ylabel) #y label
plotC
#relabel x axis:
plotD<-plotC + xlab("Pollen type") + # x axis label
  scale_x_discrete(breaks=levels(Lsmeans.df$Pollen) , labels=pollen.axis)
plotD
#Larger font:
plotE<-plotD + theme(text=element_text (size=25, face="bold") )
plotE
##remove gridlines (option), remove border
neaterA<- plotE +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.line = element_line(size = 1.5, colour = "black"))
neaterA
#axis scaling
ZeroedA<- neaterA + coord_cartesian(ylim = c(0, 5.5))
ZeroedA



#if we want y axis to start at 0:
#this removes space between x axis and bottom of bar (use for non-log scale)
JustifiedA<- neaterA +
  scale_y_continuous(expand = c(0,0), limits = c(0,5.5))
JustifiedA

fig2<-JustifiedA

