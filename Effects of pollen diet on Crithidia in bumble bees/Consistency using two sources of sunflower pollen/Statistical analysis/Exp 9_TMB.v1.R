####analysis of Exp 9_ USA_WI_SUNFOWER POLLEN EXPERIMENT ####### 
###individual bees.... JONATHAN GIACOMINI.... 4 POLLEN DIETS EFFECTS ON CRITHIDIA COUNTS#
#1 - SUN - WI Sunflower Pollen
#2 - CSUN - China Sunflower Pollen
#3 - MIX- WI Mix POllen
#4 - KMIX - Koppert Mix Pollen 

##Edited EPY 2017.06.03 : use glmmTMB to fit model and drop1() to test predictor significance

#################Content: ############################
 ##Packages
 ##Data 
 ##Models
 ##Plot
 ##Appendix: comparison with glmer.nb


####Final model -> M.nb2

##Load packages 
library(glmmTMB)
library(car)
library(bbmle) #Aic tables
library(lsmeans) #post-hoc tests with finite size correction
library(multcomp) #post-hoc tests
library(plyr)
library(dplyr)
source("glmmTMB.to.lsmeans.R") #helper function for TMB to talk to lsmeans
library(ggplot2)
library(cowplot) #multi panel plots


####Fetch data ####
##JG
D<-read.csv("/Users/rg9403/Google Drive/Sunflower Pollen Experiments/Irwin Lab_Sunflower Pollen Experiments/Exp 9_ U.S.A._WI_Sunflower Pollen Exp/U.S.A._WI_Sunflower Pollen Exp.csv",header=TRUE)
##EPY 
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/NCSU_Crithidia/Expt9")
source("glmmTMB.to.lsmeans.R") #helper function for TMB to talk to lsmeans
D<-read.csv("U.S.A._WI_Sunflower Pollen Exp.csv",header=TRUE)

str(D)
# 'data.frame':	120 obs. of  7 variables:
#   $ BeeID         : int  1 2 3 4 5 6 7 8 9 10 ...
# $ ColonyID      : Factor w/ 3 levels "A","B","C": 1 2 3 1 2 3 1 2 3 1 ...
# $ Treatment     : Factor w/ 4 levels "CSUN","KMIX",..: 1 4 3 2 1 4 3 2 1 4 ...
# $ CrithidiaCount: int  39 NA 66 127 NA 0 164 NA 6 NA ...
# $ BeeSize       : num  3.27 3.33 2.9 3.1 3 ...
# $ TimeToDeath   : int  7 1 7 7 6 7 7 5 7 1 ...
# $ Death         : int  0 1 0 0 1 0 0 1 0 1 ...


####Renaming to stay consisitent with Evan's script
Data<-dplyr::rename(D, wingmm = BeeSize, Pollen = Treatment, Count = CrithidiaCount) # renames for compatibility

#Change factor levels...
levels(Data$Pollen) ## "CSUN" "KMIX" "MIX"  "SUN" 
Data$Pollen<-factor(Data$Pollen, levels = c("KMIX", "MIX", "CSUN", "SUN"))


library(glmmTMB)
Data1<-na.omit(Data) #The following models don't work unless I omit the NAs

#Poisson:
M.p<- glmmTMB(Count ~ Pollen   + wingmm + (1|ColonyID), 
              family = "poisson",
              data=Data1) 
summary(M.p)
# AIC      BIC   logLik deviance df.resid 
# 4837.3   4853.5  -2412.7   4825.3      103 
##Deviance >> df.resid .....---> overdispersed response variable

M.nb<-glmmTMB(Count ~ Pollen   + wingmm + (1|ColonyID), 
              family = "nbinom2",
              data=Data1) 
summary(M.nb)

AICctab(M.p, M.nb) 
# dAICc  df
# M.nb    0.0 7 
# M.p  3848.2 6  
## 4K AIC units!!

#Stick with negative binomial then
drop1(M.nb, test = "Chisq")
# Df     AIC     LRT  Pr(>Chi)    
# <none>     988.79                      
# Pollen  3 1010.46 27.6760 4.248e-06 ***
#   wingmm  1  986.79  0.0048     0.945

#wingmm not helpful; recode model without it, and
  #Bring back bees without wing measures
M.nb2<-update(M.nb, ~. -wingmm, data = Data)
drop1(M.nb2, test = "Chisq")
# Df     AIC    LRT  Pr(>Chi)    
# <none>     999.08                     
# Pollen  3 1021.73 28.643 2.661e-06 ***

Lsm<-lsmeans(M.nb2, ~Pollen)
Pairs<-contrast(Lsm, "pairwise") 
Pairs
# contrast      estimate        SE df z.ratio p.value
# CSUN - KMIX -1.9377355 0.3777691 NA  -5.129  <.0001
# CSUN - MIX  -1.4659544 0.3775941 NA  -3.882  0.0006
# CSUN - SUN  -0.2328344 0.4057142 NA  -0.574  0.9399
# KMIX - MIX   0.4717811 0.3699074 NA   1.275  0.5786
# KMIX - SUN   1.7049011 0.4003308 NA   4.259  0.0001
# MIX - SUN    1.2331201 0.3942011 NA   3.128  0.0095

##Sunflowers distinguish themselves

##Get letters to plot
Lets<-cld(Lsm, Letters = letters)
Lets
Lets$.group<-gsub(" ", "", Lets$.group) ##remove spaces

################## PLOT ON SCALE OF RESPONSE ##########
#Let's exponentiate original mean+SE and mean-SE
Lsmeans.df<-as.data.frame(summary(lsmeans(M.nb2, ~Pollen)))
Lsmeans.df$Count<-exp(Lsmeans.df$lsmean)
Lsmeans.df$btlo<-exp(Lsmeans.df$lsmean - Lsmeans.df$SE) #lower bound
Lsmeans.df$bthi<-exp(Lsmeans.df$lsmean + Lsmeans.df$SE) #upper bound

Lsmeans.df

#Merge with letters
Lets2<- dplyr::select(Lets, Pollen, .group)
Lets2

Lsmeans.df.merged<- join(Lsmeans.df, Lets2, by = "Pollen")

Lsmeans.df.merged
# Pollen   lsmean        SE df asymp.LCL asymp.UCL    Count      btlo      bthi .group
# 1   CSUN 2.563754 0.4287006 NA  1.723517  3.403992 12.98448  8.457502  19.93456      a
# 2   KMIX 4.501490 0.4197674 NA  3.678761  5.324219 90.15135 59.247437 137.17499      b
# 3    MIX 4.029709 0.4209153 NA  3.204730  4.854688 56.24453 36.921482  85.68041      b
# 4    SUN 2.796589 0.4465523 NA  1.921362  3.671815 16.38865 10.485953  25.61405      a

Lsmeans.df<-Lsmeans.df.merged #overwrites old

#x-axis:
pollen.axis<-c("WF Mix 1","WF Mix 2","Sun (CN)", "Sun (US)")
#y-axis:
ylabel<- expression(paste("Parasite load (cells * 0.02 ", mu, L^-1, "))", sep="") ) #check for encoding

#ready to plot?
p<-ggplot(Lsmeans.df, aes(x = Pollen, y = Count)) + 
  #Option to sprinkle raw data... 
  geom_point(data = Data, aes( x = Pollen, y = Count), position = "jitter", 
             color = "blue", shape = 25, size = 0.5, alpha = 0.4)+
  geom_pointrange(aes(ymin = btlo, ymax = bthi), size = 1, color = "black") +
  scale_x_discrete(labels=pollen.axis) +
  ylab(ylabel) + #y label
  xlab("Pollen")
p
#Add text
Panel.1c <- p + geom_text(aes(y = bthi + 0.2 * max(bthi), label = .group),
                          size = 8)
Panel.1c
# ggsave("Exp9_v2.TMB.pdf", height = 5, width =5)


################# Composite figures #####################
#Use  plot_grid

#Must run other scripts or load workspaces:
#setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Analysis.Individual.Bee.and.Farm/Analyses/Individual Bees/Crithidia Load")
#source("negative_binomial_individuals_v2.TMB.R")
#setwd("C:/Users/Evan/Google Drive/Project_Sunflower/NCSU_Crithidia/Exp2_JCRoulstonStrain_BuckSunKMix_2Week")
#source("Exp2.buck.kopp.sun.TMB.v1a.R")


Panel.ABC<-plot_grid(Panel.1a, Panel.1b, Panel.1c)
Panel.ABC

##Notes: Panels are 
#(A) Monofloral vs Mixed pollen
#(B) Effect of diet post-infection
#(C) Consistency of sunflower effect 
#(Supp) Reproducibility of effect... former panel C

#To-do:
##Add Panel Titles
##Remove x axis title
##Remove redundant y-axis titles

F1a.titled <- 
  Panel.1a + 
  ggtitle("(A) Monofloral vs mixed pollen")+
  theme(axis.title.x = element_blank())
F1b.titled <- 
  Panel.1b + 
  ggtitle("(B) Effect of diet post-infection")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
F1c.titled <- 
  Panel.1c + 
  ggtitle("(C) Consistency of sunflower effect")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

F1.all<- plot_grid(F1a.titled, F1b.titled, F1c.titled,
                   nrow = 1, rel_widths = c(1.2, 1, 1.2), align = "h")
F1.all

#where to put this??
setwd("C:/Users/Evan/Google Drive/Project_Sunflower")
ggsave(F1.all, file = "Fig.1.3panel.pdf", height = 4, width = 12)



###################### End main script #############################################

###########################################################################################

##################### Appendix : Compare with lme4 glmer.nb fit ######################


###try with a different package for comparison
library(lme4)
M0<- glmer.nb(totcrith ~ Pollen  + wingmm + (1|ColonyID)  , data=Data)

summary(M0) #Same thing...wingmm not significnat p = 0.945690

Anova(M0)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: totcrith
# Chisq Df Pr(>Chisq)    
# Pollen 34.2150  3  1.785e-07 ***
#   wingmm  0.0046  1     0.9457

# Let's run it again without wingmm
M2<- glmer.nb(totcrith ~ Pollen  + (1|ColonyID)  , data=Data)

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

##Similar to TMB model (note different order of terms :))
Lsm<-lsmeans(M.nb2, ~Pollen)
Pairs<-contrast(Lsm, "pairwise") 
Pairs
# contrast      estimate        SE df z.ratio p.value
# KMIX - MIX   0.4717837 0.3699073 NA   1.275  0.5786
# KMIX - CSUN  1.9377381 0.3777690 NA   5.129  <.0001
# KMIX - SUN   1.7049027 0.4003308 NA   4.259  0.0001
# MIX - CSUN   1.4659544 0.3775938 NA   3.882  0.0006
# MIX - SUN    1.2331189 0.3942008 NA   3.128  0.0095
# CSUN - SUN  -0.2328354 0.4057139 NA  -0.574  0.9399


