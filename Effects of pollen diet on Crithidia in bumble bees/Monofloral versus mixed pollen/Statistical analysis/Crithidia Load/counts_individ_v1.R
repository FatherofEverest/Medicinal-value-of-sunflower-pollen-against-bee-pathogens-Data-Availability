####analysis of POLLEN EXPERIMENT ####### 
###individual bees.... JONATHAN GIACOMINI.... 4 POLLEN DIETS EFFECTS ON CRITHIDIA COUNTS#


####Final model -> nbinom0 

#Jon-- I will try to guide you through this
rm(list=ls()) #clear memory
#set working directory
setwd() #choose source file location, if desired

#read in Jon's data from somewhere
#data<-read.csv("/Users/rg9403/Desktop/PollenEXP/PollenExp_Winter2014/Data_CSVs/Analysis_1_Crithidia load.csv",header=TRUE)
#epy machine:
data<-read.csv('C:/Users/Evan/Google Drive/ The Sunflower Pollen Project/DATA_ANALYSES Individual Bee +  Farm/DATA/CSV/Analysis_1_Crithidia load.csv')
str(data)

#### Check correlations between explanatory variables, especially collectiondate.####
d<- data.frame(data$DissectWeight,data$WingSize,data$InocDate,data$CallowMass)
cor(d,use="complete.obs")
####WingSize and DissectWeight and callowmass 
#are highly correlated, 
#thus we should include at most ONE OF THE THREE in model


####Renaming to stay consisitent with Evan's script
Data<-data #to rename 
#Change InocDate to factor
Data$InocDate<-factor(Data$InocDate)
Data$wingmm<-Data$WingSize # renames for compatibility
Data$Pollen<-Data$Treatment #ditto
Data$totcrith<-Data$Count #ditto

#install.packages("R2admb")
#install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                         getOption("repos")), type="source")
library(glmmADMB)
View(Data)
length(data$BeeID)
Data1<-na.omit(Data) #The following models don't work unless I omit the NAs
length(Data1$BeeID) #removed 2 bees, no big deal

#scout distribution, see helpful script by julia pilowsky
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
#that's quite OK

Boxplot(Data1$Count~Data1$Treatment)
#variance increasing with mean
#negative binomial should account for that

#Compare Poisson with zero inflation, regular negative binomial, and negative binomial with zero inflation
#first try poisson model:
fitPoiss<- glmmadmb(totcrith ~ Pollen  + wingmm + (1|ColonyID)  + (1|InocDate), 
                    family = "poisson",
                    data=Data1)  
summary(fitPoiss)

#######################################  NEG BIN ########################################################
#####now negative binomial
#Start with full nesting

nbinom0<-glmmadmb(totcrith ~ Pollen  + wingmm +(1|ColonyID)  + (1|InocDate), 
                  family = "nbinom",
                  data=Data1)  
library(bbmle)
AICtab(fitPoiss, nbinom0)
summary(nbinom0) #get some coefficients 
# Call:
#   glmmadmb(formula = totcrith ~ Pollen + wingmm + (1 | ColonyID) + 
#              (1 | InocDate), data = Data1, family = "nbinom")
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

#deviance table (chi squared tests)
Anova(nbinom0)

library(lsmeans)
lsmeans(nbinom0, ~Pollen)

library(multcomp)
#post hoc test 
summary(glht(nbinom0,linfct=mcp(Pollen="Tukey")))
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Tukey Contrasts
# 
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
# 2     8 -874.68  1    6.644 0.009949 ** 
### We can't drop Wingsize fromt he model

#compare AIC values (lower is better)
#more than 2 aic units is usually called "significant difference"
library(bbmle)
AIC(nbinom0) 
#[1] 1765.37
AIC(m1) #this does the same, but manually
#[1] 1770.014

#### WingSize and DissectWeight are highly correlated, 
#so let's swap them in the model and see what is a better fit
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
  geom_bar(position=position_dodge(), 
           stat="identity", colour="black", size=1.5 ) +
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE),
                size=1.5, width=0.2,  # Width of the error bars
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
#find my axes!
axed<-neater+ theme(axis.line.y=element_line(size=2, color="black"),
              axis.line.x=element_line(size=2, color="black"))
axed + coord_cartesian(ylim = c(0, 6.5))
#prefer "axed"

#if we want y axis to start at 0:
#this removes space between x axis and bottom of bar (use for non-log scale)
#Justified<- neater +
 # scale_y_continuous(expand = c(0,0), limits = c(0,6.5))
#Justified
#adjust legend
#fig1<-MoveLeg #or "Justified" to start y axis at 0
#fig1<-fig1 +  theme(axis.title.y=element_text(vjust=0.9))
fig1<-axed
#export to PDF
#dev.new() #sometimes this helps with errors
setwd("C:/Users/Evan/Google Drive/ The Sunflower Pollen Project/DATA_ANALYSES Individual Bee +  Farm/Analyses")
pdf("fig1_counts_v2.pdf",height=7, width=7, paper='special') 
fig1
dev.off()

jpeg("fig1_counts_v2.jpg")
fig1
dev.off()

#for eps
setEPS()
postscript("fig1countsv2.eps", horizontal = FALSE, onefile = FALSE, 
           paper = "special", height = 7, width = 7)
fig1
dev.off()





#BUCK
exp(4.6954957) #109.4531 x50 = 5472.655 cells/uL
#LSmeans +/- LS_SEs
4.6954957+1.688676 #= 6.384172
4.6954957-1.688676 #= 3.00682
#Backtransformed SEs
exp(6.384172) #592.394 * 50 = 29619.7 cells/uL
exp (3.00682) #20.22299 * 50 = 1011.149 cells/uL



#SUN
exp(0.7166957) #2.047656 x 50 = 102.3828 cells/uL
#LSmeans +/- LS_SEs
0.7166957+1.756342 #2.473038
0.7166957-1.756342 #-1.039646
#Backtransformed SEs
exp(2.473038) #11.85842 * 50 = 592.921 cells/uL
exp(-1.039646) #0.3535798 * 50 = 17.67899 cells/uL

#RAPE
exp(3.7538951) #42.68703 *50 = 2134.352 cells/uL
#LSmeans +/- LS_SEs
3.7538951+1.701106 #= 5.455001
3.7538951-1.701106 #= 2.052789
#Backtransformed SEs
exp(5.455001) #233.9251 * 50 = 11696.25 cells/uL
exp (2.052789) #7.789596 * 50 = 389.4798 cells/uL

#MIX
exp(2.5084296) #12.28562 *50 = 614.281 cells/uL
#LSmeans +/- LS_SEs
2.5084296+1.744933 #= 4.253363
2.5084296-1.744933 #= 0.763496
#Backtransformed SEs
exp(4.253363) #70.34157 * 50 = 3517.079 cells/uL
exp (0.763496) #2.145765 * 50 = 107.2882 cells/uL