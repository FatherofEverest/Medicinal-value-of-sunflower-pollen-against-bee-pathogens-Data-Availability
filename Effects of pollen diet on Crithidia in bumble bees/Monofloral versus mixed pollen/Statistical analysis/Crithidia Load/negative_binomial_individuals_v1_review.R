####analysis of POLLEN EXPERIMENT ####### 
###individual bees.... JONATHAN GIACOMINI.... 4 POLLEN DIETS EFFECTS ON CRITHIDIA COUNTS#

library(cowplot)
library(ggplot2)
library(glmmADMB)
library(lme4)
library(plyr)
library(dplyr)
library(lsmeans)

####Final model -> nbinom0 

#set working directory
rm(list=ls()) #clear memory
setwd() #choose source file location, if desired

#EPY:
#setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Analysis.Individual.Bee.and.Farm/Analyses/Individual Bees/Crithidia Load")

#read in Jon's data from somewhere
#EPY:
data<-read.csv("Analysis_1_Crithidia load.csv",header=TRUE)
str(data)
data<-read.csv("/Users/rg9403/Desktop/PollenEXP/PollenExp_Winter2014/Data_CSVs/Analysis_1_Crithidia load.csv", header=TRUE)
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

#install.packages("R2admb")
# install.packages("glmmADMB", 
#                  repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                          getOption("repos")),
#                  type="source")
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

#Compare the 2 models:
anova(nbinom0, zinfl)
# Model 1: totcrith ~ Pollen + wingmm
# Model 2: totcrith ~ Pollen + wingmm
# NoPar  LogLik Df Deviance Pr(>Chi)
# 1     8 -874.68                     
# 2     9 -874.68  1    0.002   0.9643
  #zero-inflation not helpful

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
library(lsmeans)
library(ggplot2)
#Use lsmeans to get mean counts
My.Lsmeans <- lsmeans(nbinom0, ~Pollen) 
My.Lsmeans
cld(My.Lsmeans, Letters = LETTERS)
#This will give significance dislays that you can add to the plot with geom_text
# Pollen    lsmean       SE df  asymp.LCL asymp.UCL .group
# Sun    0.7166957 1.756342 NA -2.7256715  4.159063  A    
# Mix    2.5084296 1.744933 NA -0.9115772  5.928436   B   
# Rape   3.7538951 1.701106 NA  0.4197876  7.088003    C  
# Buck   4.6954957 1.688676 NA  1.3857523  8.005239     D 
# 

#Run pairwise contrasts
Contrasts<- lsmeans::contrast(My.Lsmeans, "pairwise")
Contrasts
#make a data frame from the lsmeans object
Lsmeans.df<-summary(My.Lsmeans)
str(Lsmeans.df)
#View(Lsmeans.df)
#label for x axis
pollen.axis<-c( "Buckwheat","Mix","Rape", "Sunflower")
plot0<- ggplot(Lsmeans.df, aes(x=Pollen, y=lsmean)) + 
  geom_bar(position=position_dodge(), stat="identity", colour="black" ) +
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE),
                size=1, width=0.2,  # Width of the error bars
                position=position_dodge(0.9))
plot0 #Solid base plot
plot01<-plot0 + theme_bw() #remove colored background
ylabel<- expression(bold(paste("Parasite load (ln (cells * 0.02 ", mu, L^-1, "))", sep="") )) #check for encoding
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


#if we want y axis to start at 0:
#this removes space between x axis and bottom of bar (use for non-log scale)
# Justified<- neater +
#   scale_y_continuous(expand = c(0,0), limits = c(0,6.5))
# Justified
#adjust legend
#fig1<-MoveLeg #or "Justified" to start y axis at 0
#fig1<-fig1 +  theme(axis.title.y=element_text(vjust=0.9))
fig1<-neater
#export to PDF
#dev.new() #sometimes this helps with errors
# setwd("C:/Users/Evan/Dropbox/Jess Sum15/jess stats/Counts/R_analyses")
# pdf("fig1_counts_v2.pdf",height=7, width=7, paper='special') 
# fig1
# dev.off()

#for eps
# setEPS()
# postscript("fig1countsv2.eps", horizontal = FALSE, onefile = FALSE, 
#            paper = "special", height = 7, width = 7)
# fig1
# dev.off()


################## OPTION 2: PLOT ON SCALE OF RESPONSE ##########
#Let's exponentiate original mean+SE and mean-SE
Lsmeans.df<-as.data.frame(summary(lsmeans(nbinom0, ~Pollen)))
Lsmeans.df$btmean<-exp(Lsmeans.df$lsmean)
Lsmeans.df$btlo<-exp(Lsmeans.df$lsmean - Lsmeans.df$SE) #lower bound
Lsmeans.df$bthi<-exp(Lsmeans.df$lsmean + Lsmeans.df$SE) #upper bound

#rename "btmean" to match name of counts in original data...
 #will be important when adding points to the plot
names(Data1)
library(plyr)
library(dplyr)
Lsmeans.df<-dplyr::rename(Lsmeans.df, Count = lsmean)
Lsmeans.df

#compare to asking for response scale initially:
lsmeans(nbinom0, ~Pollen, type = "response")
#ready to plot?
p<-ggplot(Lsmeans.df, aes(x = Pollen, y = Count)) + 
  geom_pointrange(aes(ymin = btlo, ymax = bthi))
p
  

#Option to sprinkle raw data... 
p + geom_point(data = Data1, aes( x = Pollen, y = Count), color = "blue", shape = 25)
#error bars seem quite big... not sure what's happening

confint(nbinom0) #gives 2.5% and 97.5% bounds
# (Intercept)  5.309049 10.0760729
# PollenMix   -2.802730 -1.5714022
# PollenRape  -1.531138 -0.3520631
# PollenSun   -4.651322 -3.3062775
# wingmm      -3.190393 -0.4450718
Lsmeans.df
# Pollen     Count       SE df  asymp.LCL asymp.UCL     btmean       btlo      bthi
# 1   Buck 4.6954957 1.688676 NA  1.3857523  8.005239 109.453052 20.2229904 592.39363
# 2    Mix 2.5084296 1.744933 NA -0.9115772  5.928436  12.285622  2.1457649  70.34158
# 3   Rape 3.7538951 1.701106 NA  0.4197876  7.088003  42.687028  7.7895928 233.92524
# 4    Sun 0.7166957 1.756342 NA -2.7256715  4.159063   2.047656  0.3535797  11.85842

#It looks like a lot of the uncertainty is coming from the "intercept" piece
#where CI spans 5.3 - 10.0 ON LOG SCALE!!!

#If one wants to plot the fixed effects confidence interval instead:
CI.df<-confint(nbinom0)
Fixed<-fixef(nbinom0)
str(CI.df)
CI.and.est<-as.data.frame(cbind(CI.df, Fixed))
CI.and.est
CI.and.est$Effect<-rownames(CI.and.est)
CI.and.est
names(CI.and.est)
ggplot(CI.and.est, aes(x=Effect, y = Fixed, ymin = `2.5 %`, ymax = `97.5 %`))+
  geom_pointrange()+ ylab("Estimate") + xlab ("Term")

#axis scaling; don't use for log scale
# Zeroed<- neater + coord_cartesian(ylim = c(0, 6.5))
# Zeroed


#Try to diagnose again:
m.nb<-glmmadmb(totcrith ~ 
                         Pollen  + wingmm +(1|ColonyID)  + (1|InocDate), 
                       family = "nbinom",
                       data=Data1)  #got an error with Jess's data
summary(m.nb)
lsm2<-lsmeans(m.nb, ~Pollen, data = Data1)
plot(lsm2) #still quite spread out, doesn't match fixed effects

pairs.2<-contrast(lsm2, method ="pairwise")
pairs.2




#citation info
# setwd("C:/Users/Evan/Dropbox/Jess Sum15/jess stats/Counts/R_analyses")
# #install.packages("bibtex")
# library(bibtex)
# write.bib(entry="ggplot2", file = "ggplot2.bib", append = FALSE, verbose = TRUE)
# write.bib(entry="coxme", file = "coxme.bib", append = FALSE, verbose = TRUE)
# write.bib(entry="survival", file = "survival.bib", append = FALSE, verbose = TRUE)
# write.bib(entry="car", file = "car.bib", append = FALSE, verbose = TRUE)
# write.bib(entry="lme4", file = "lme4.bib", append = FALSE, verbose = TRUE)
# write.bib(entry="MASS", file = "MASS.bib", append = FALSE, verbose = TRUE)
# write.bib(entry="lsmeans", file = "lsmeans.bib", append = FALSE, verbose = TRUE)
# write.bib(entry="glmmADMB", file = "glmmADMB.bib", append = FALSE, verbose = TRUE)
# 





