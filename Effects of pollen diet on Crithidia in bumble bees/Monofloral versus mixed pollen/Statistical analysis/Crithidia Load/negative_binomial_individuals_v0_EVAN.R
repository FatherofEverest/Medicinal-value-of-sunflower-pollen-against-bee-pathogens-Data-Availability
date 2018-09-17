####analysis of POLLEN EXPERIMENT ####### 
###individual bees.... JONATHAN GIACOMINI.... 4 POLLEN DIETS EFFECTS ON CRITHIDIA COUNTS#


####Final model -> nbinom0 

#Jon-- I will try to guide you through this
rm(list=ls()) #clear memory
#set working directory
setwd() #choose source file location, if desired
#read in Jon's data from somewhere
data<-read.csv("/Users/rg9403/Desktop/PollenEXP/PollenExp_Winter2014/Data_CSVs/Analysis_1_Crithidia load.csv",header=TRUE)
str(data)
# 'data.frame':	234 obs. of  8 variables:
#   $ BeeID        : int  1 4 5 6 7 8 13 15 17 18 ...
# $ ColonyID     : Factor w/ 8 levels "JG 1","JG 2",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Treatment    : Factor w/ 4 levels "Buck","Mix","Rape",..: 4 2 4 1 3 2 4 3 4 1 ...
# $ Count        : int  14 1 1 5 85 60 2 500 1 109 ...
# $ DissectWeight: num  0.15 0.123 0.0905 0.137 0.1035 ...
# $ WingSize     : num  1.76 1.76 1.44 1.76 1.48 1.36 1.44 1.28 1.2 1.52 ...
# $ InocDate     : int  1 1 1 2 2 2 3 3 6 9 ...
# $ CallowMass   : num  0.15 0.151 0.097 0.161 0.093 ...

#### Check correlations between explanatory variables, especially collectiondate.####
d<- data.frame(data$DissectWeight,data$WingSize,data$InocDate,data$CallowMass)
cor(d,use="complete.obs")

#                         data.DissectWeight data.WingSize data.InocDate data.CallowMass
# data.DissectWeight        1.000000000    0.87582779   0.002730283     0.879884830
# data.WingSize             0.875827794    1.00000000  -0.056381693     0.949437393
# data.InocDate             0.002730283   -0.05638169   1.000000000    -0.003681552
# data.CallowMass           0.879884830    0.94943739  -0.003681552     1.000000000
####WingSize and DissectWeight are highly correlated, thus we should toss one of them from the model


####Renaming to stay consisitent with Evan's script
Data<-data #to rename 
#Change InocDate to factor
Data$InocDate<-factor(Data$InocDate)
Data$wingmm<-Data$WingSize # renames for compatibility
Data$Pollen<-Data$Treatment #ditto
Data$totcrith<-Data$Count #ditto

install.packages("R2admb")
install.packages("glmmADMB", 
                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
                         getOption("repos")),
                 type="source")
library(glmmADMB)
Data1<-na.omit(Data) #The following models don't work unless I omit the NAs
#Compare Poisson with zero inflation, regular negative binomial, and negative binomial with zero inflation
#Example:
#fit_zipoiss <- glmmadmb(NCalls ~ (FoodTreatment + ArrivalTime) * SexParent +
#offset(log(BroodSize)) + (1 | Nest), data = Owls, zeroInflation = TRUE,
#family = "poisson") 
#first try poisson model:
fitPoiss<- glmmadmb(totcrith ~ Pollen  + wingmm + (1|ColonyID)  + (1|InocDate), 
                    family = "poisson",
                    data=Data1)  #my data set was throwing errors with poisson fit, no good
summary(fitPoiss)
# Call:
#   glmmadmb(formula = totcrith ~ Pollen + wingmm + (1 | ColonyID) + 
#              (1 | InocDate), data = Data1, family = "poisson")
# 
# AIC: 4168.2 
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  11.3302     0.5971   18.98   <2e-16 ***
#   PollenMix    -0.7493     0.0841   -8.91   <2e-16 ***
#   PollenRape    0.1853     0.0823    2.25    0.024 *  
#   PollenSun    -2.7861     0.1407  -19.80   <2e-16 ***
#   wingmm       -5.7518     0.3007  -19.13   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Number of observations: total=232, ColonyID=8, InocDate=37 
# Random effect variance(s):
#   Group=ColonyID
# Variance StdDev
# (Intercept)   0.5047 0.7104
# Group=InocDate
# Variance StdDev
# (Intercept)    1.957  1.399
# 
# 
# Log-likelihood: -2077.1 

#####try zero inflation with negative binomial:
zinfl<- glmmadmb(totcrith ~ Pollen  + wingmm + (1|ColonyID)  + (1|InocDate), 
                                    family = "nbinom",zeroInflation = TRUE, data=Data1) 
# got error with jess's data 
summary(zinfl)
# Call:
#   glmmadmb(formula = totcrith ~ Pollen + wingmm + (1 | ColonyID) + 
#              (1 | InocDate), data = Data1, family = "nbinom", zeroInflation = TRUE)
# 
# AIC: 1767.4 
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    7.692      1.216    6.33  2.5e-10 ***
#   PollenMix     -2.187      0.314   -6.96  3.3e-12 ***
#   PollenRape    -0.942      0.301   -3.13   0.0017 ** 
#   PollenSun     -3.979      0.343  -11.60  < 2e-16 ***
#   wingmm        -1.818      0.700   -2.60   0.0095 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Number of observations: total=232, ColonyID=8, InocDate=37 
# Random effect variance(s):
#   Group=ColonyID
# Variance StdDev
# (Intercept)   0.1898 0.4357
# Group=InocDate
# Variance StdDev
# (Intercept)   0.2708 0.5204
# 
# Negative binomial dispersion parameter: 0.45347 (std. err.: 0.050877)
# Zero-inflation: 1.0641e-06  (std. err.:  0.00019545 )
# 
# Log-likelihood: -874.684 


#######################################  NEG BIN ########################################################
#####now negative binomial
#Start with full nesting

nbinom0<-glmmadmb(totcrith ~ Pollen  + wingmm +(1|ColonyID)  + (1|InocDate), 
                  family = "nbinom",
                  data=Data1)  #got an error with Jess's data
summary(nbinom0) #get some coefficients 
# Call:
#   glmmadmb(formula = totcrith ~ Pollen + wingmm + (1 | ColonyID) + 
#              (1 | InocDate), data = Data1, family = "nbinom")
# 
# AIC: 1765.4 
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    7.693      1.216    6.33  2.5e-10 ***
#   PollenMix     -2.187      0.314   -6.96  3.3e-12 ***
#   PollenRape    -0.942      0.301   -3.13   0.0017 ** 
#   PollenSun     -3.979      0.343  -11.60  < 2e-16 ***
#   wingmm        -1.818      0.700   -2.60   0.0094 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Number of observations: total=232, ColonyID=8, InocDate=37 
# Random effect variance(s):
#   Group=ColonyID
# Variance StdDev
# (Intercept)   0.1899 0.4358
# Group=InocDate
# Variance StdDev
# (Intercept)   0.2708 0.5203
# 
# Negative binomial dispersion parameter: 0.45344 (std. err.: 0.050853)
# 
# Log-likelihood: -874.685 

library(multcomp)
#post hoc test 
summary(glht(nbinom0,linfct=mcp(Pollen="Tukey")))
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Tukey Contrasts
# 
# 
# Fit: glmmadmb(formula = totcrith ~ Pollen + wingmm + (1 | ColonyID) + 
#                 (1 | InocDate), data = Data1, family = "nbinom")
# 
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# Mix - Buck == 0   -2.1871     0.3141  -6.963  < 0.001 ***
#   Rape - Buck == 0  -0.9416     0.3008  -3.130  0.00845 ** 
#   Sun - Buck == 0   -3.9788     0.3431 -11.596  < 0.001 ***
#   Rape - Mix == 0    1.2455     0.4243   2.935  0.01508 *  
#   Sun - Mix == 0    -1.7917     0.4920  -3.641  0.00125 ** 
#   Sun - Rape == 0   -3.0372     0.4953  -6.132  < 0.001 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)


#you can do a shortcut to test individual terms like this:
m1<-update(nbinom0, ~. - wingmm) #this reminds me of how they take taco bell orders,
    # [my taco] -chz -sauce
summary(m1)
anova(nbinom0, m1) #if p>0.05, you are good to drop the parameter (in this case wingmm)
                 ####and keep the simpler model
# Analysis of Deviance Table
# 
# Model 1: totcrith ~ Pollen
# Model 2: totcrith ~ Pollen + wingmm
# NoPar  LogLik Df Deviance Pr(>Chi)   
# 1     7 -878.01                        
# 2     8 -874.68  1    6.644 0.009949 ** ### We can't drop Wingsize fromt he model

library(car)
Anova(m1) #this will give chisquare tests, similar to an anova table, called "Deviance Table"

#compare AIC values (lower is better)
#more than 2 aic units is usually called "significant difference"
library(bbmle)
AIC(nbinom0) 
#[1] 1765.37
AIC(m1) #this does the same, but manually
#[1] 1770.014

#### WingSize and DissectWeight are highly correlated, so let's swap them in the model and see what is a better fit
m2<-update(nbinom0, ~. - wingmm + DissectWeight)
summary(m2)

m3<-update(nbinom0, ~. - wingmm + CallowMass)
summary(m3)

AIC(m2)
# [1] 1768.058
AIC(m1)
#[1] 1770.014
AIC(nbinom0) 
#[1] 1765.37  #### WINNER WINNER, CHICKEN DINNER!
AIC(m3)
# [1] 1765.566

#2. likelihood ratio test (alternative to chi-squared)
Stripped<-update(nbinom0, ~. - wingmm) #makes a reduced model
anova(nbinom0, Stripped) #compares reduced model to the model with pollen as predictor
# Analysis of Deviance Table
# 
# Model 1: totcrith ~ Pollen
# Model 2: totcrith ~ Pollen + wingmm
# NoPar  LogLik Df Deviance Pr(>Chi)   
# 1     7 -878.01                        
# 2     8 -874.68  1    6.644 0.009949 **




###try with a different package for comparison
library(lme4)
M0<- glmer.nb(totcrith ~ Pollen  + wingmm + (1|ColonyID)  + (1|InocDate), data=Data1)
#warning message-- model not converging? 
summary(M0)
Anova(M0)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: totcrith
# Chisq Df Pr(>Chisq)    
# Pollen 425443  3  < 2.2e-16 ***
#   wingmm  65967  1  < 2.2e-16 ***

#you can do the same model simplification workflow
m1a<-update(M0, ~. - wingmm) #this reminds me of how they take taco bell orders,
# [my taco] -chz -sauce
summary(m1a)
anova(M0, m1a) #if p>0.05, you are good to drop the parameter (in this case wingmm)
####and keep the simpler model
summary(m1a)
Anova(m1a) #this will give chisquare tests, similar to an anova table, called "Deviance Table"
Stripped<-update(m1a, ~. - Pollen) #makes a reduced model
anova(m1a, Stripped) #compares reduced model to the model with pollen as predictor


###################################### PLOTS ####################################################
#Use lsmeans to get mean counts
library(lsmeans)
My.Lsmeans <- lsmeans(nbinom0, ~Pollen) 
My.Lsmeans
altLsmeans<- lsmeans(M0, ~Pollen, data=Data1) #substitute the name of your favorite model
#choose the simplest model you can without excluding significant terms
#obviously error messages may narrow down your choices:)
altLsmeans


###HERE IS SOME UNEDITED SCRIPT TO PLAY WITH:
####PLOTTING
#Use lsmeans to get mean counts
library(lsmeans)
My.Lsmeans <- lsmeans(CountFinal, ~Pollen)
My.Lsmeans
altLsmeans<- lsmeans(M1, ~Pollen)
altLsmeans #somewhat smaller confidence bands
#I am going to go with the glmmadmb model ("CountFinal"), 
#since it did not give any warnings

#BELOW adapted FROM CONROY EXPERIMENT
#######PLOTTING : PARASITE LOAD#############
#plotting marginal means:
library(ggplot2)
#make a data frame from the lsmeans object
Lsmeans.df<-summary(My.Lsmeans)
str(Lsmeans.df)
View(Lsmeans.df)
#label for x axis
pollen.axis<-c( "Buckwheat","Mix","Rape", "Sunflower")

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
Zeroed<- neater + coord_cartesian(ylim = c(0, 6.5))
Zeroed



#if we want y axis to start at 0:
#this removes space between x axis and bottom of bar (use for non-log scale)
Justified<- neater +
  scale_y_continuous(expand = c(0,0), limits = c(0,6.5))
Justified
#adjust legend
#fig1<-MoveLeg #or "Justified" to start y axis at 0
#fig1<-fig1 +  theme(axis.title.y=element_text(vjust=0.9))
fig1<-Justified
#export to PDF
#dev.new() #sometimes this helps with errors
setwd("C:/Users/Evan/Dropbox/Jess Sum15/jess stats/Counts/R_analyses")
pdf("fig1_counts_v2.pdf",height=7, width=7, paper='special') 
fig1
dev.off()

#for eps
setEPS()
postscript("fig1countsv2.eps", horizontal = FALSE, onefile = FALSE, 
           paper = "special", height = 7, width = 7)
fig1
dev.off()

#citation info
setwd("C:/Users/Evan/Dropbox/Jess Sum15/jess stats/Counts/R_analyses")
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







#### TREATMENT ####
library(lme4)
mod1<-glmer(Count ~ Treatment+(1|BeeID)+(1|ColonyID),family=poisson,data=data)
summary(mod1)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: poisson  ( log )
# Formula: Count ~ Treatment + (1 | BeeID) + (1 | ColonyID)
# Data: data
# 
# AIC      BIC   logLik deviance df.resid 
# 1794.8   1815.5   -891.4   1782.8      228 
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -1.07805 -0.43158 -0.00249  0.04370  0.38878 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# BeeID    (Intercept) 3.4560   1.859   
# ColonyID (Intercept) 0.2062   0.454   
# Number of obs: 234, groups:  BeeID, 234; ColonyID, 8
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)     3.6461     0.3020  12.072  < 2e-16 ***
#   TreatmentMix   -2.5168     0.3740  -6.729  1.7e-11 ***
#   TreatmentRape  -0.7533     0.3573  -2.109    0.035 *  
#   TreatmentSun   -4.4499     0.3962 -11.232  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) TrtmnM TrtmnR
# TreatmentMx -0.539              
# TreatmentRp -0.571  0.458       
# TreatmentSn -0.497  0.433  0.431


#### INOCULATION DATE #####
mod2<-glmer(Count ~ InocDate+(1|BeeID)+(1|ColonyID),family=poisson,data=data)
summary(mod2)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: poisson  ( log )
# Formula: Count ~ InocDate + (1 | BeeID) + (1 | ColonyID)
# Data: data
# 
# AIC      BIC   logLik deviance df.resid 
# 1920.4   1934.3   -956.2   1912.4      230 
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -0.66219 -0.59726  0.02156  0.04200  0.05782 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# BeeID    (Intercept) 6.640    2.577   
# ColonyID (Intercept) 0.117    0.342   
# Number of obs: 234, groups:  BeeID, 234; ColonyID, 8
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  1.96710    0.51975   3.785 0.000154 ***
#   InocDate    -0.01480    0.01783  -0.830 0.406439    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr)
# InocDate -0.887


#### WING SIZE ####
mod3<-glmer(Count ~ WingSize+(1|BeeID)+(1|ColonyID),family=poisson,data=data)
summary(mod3)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: poisson  ( log )
# Formula: Count ~ WingSize + (1 | BeeID) + (1 | ColonyID)
# Data: data
# 
# AIC      BIC   logLik deviance df.resid 
# 1914.2   1928.1   -953.1   1906.2      229 
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -0.70036 -0.57254  0.01137  0.04127  0.06779 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# BeeID    (Intercept) 6.5019   2.5499  
# ColonyID (Intercept) 0.1614   0.4017  
# Number of obs: 233, groups:  BeeID, 233; ColonyID, 8
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)    5.154      1.757   2.933  0.00336 **
#   WingSize      -2.152      1.059  -2.033  0.04206 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr)
# WingSize -0.990


#### TREATMENT + WING SIZE ####
mod4<-glmer(Count ~ Treatment+WingSize+(1|BeeID)+(1|ColonyID),family=poisson,data=data)
summary(mod4)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: poisson  ( log )
# Formula: Count ~ Treatment + WingSize + (1 | BeeID) + (1 | ColonyID)
# Data: data
# 
# AIC      BIC   logLik deviance df.resid 
# 1784.7   1808.8   -885.3   1770.7      226 
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -1.05944 -0.40502 -0.01219  0.04313  0.46337 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# BeeID    (Intercept) 3.2807   1.8113  
# ColonyID (Intercept) 0.2752   0.5246  
# Number of obs: 233, groups:  BeeID, 233; ColonyID, 8
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)     7.7986     1.3624   5.724 1.04e-08 ***
#   TreatmentMix   -2.5238     0.3673  -6.870 6.40e-12 ***
#   TreatmentRape  -0.8736     0.3507  -2.491  0.01274 *  
#   TreatmentSun   -4.5197     0.3902 -11.582  < 2e-16 ***
#   WingSize       -2.4883     0.7968  -3.123  0.00179 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) TrtmnM TrtmnR TrtmnS
# TreatmentMx -0.175                     
# TreatmentRp -0.231  0.459              
# TreatmentSn -0.215  0.433  0.437       
# WingSize    -0.974  0.061  0.111  0.111

library(car)
Anova(mod4)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: Count
#             Chisq    Df   Pr(>Chisq)    
# Treatment   151.3553  3   < 2.2e-16 ***
# WingSize    8.7123    1   0.003161 ** 

library(multcomp)
#post hoc test 
summary(glht(mod4,linfct=mcp(Treatment="Tukey")))
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Tukey Contrasts
# 
# 
# Fit: glmer(formula = Count ~ Treatment + WingSize + (1 | BeeID) + 
#              (1 | ColonyID), data = data, family = poisson)
# 
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# Mix - Buck == 0   -2.5238     0.3673  -6.870   <0.001 ***
#   Rape - Buck == 0  -0.8736     0.3507  -2.491   0.0612 .  
# Sun - Buck == 0   -4.5197     0.3902 -11.582   <0.001 ***
#   Rape - Mix == 0    1.6502     0.3736   4.417   <0.001 ***
#   Sun - Mix == 0    -1.9959     0.4039  -4.942   <0.001 ***
#   Sun - Rape == 0   -3.6461     0.3946  -9.240   <0.001 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)

plot(data$WingSize,data$Count)
