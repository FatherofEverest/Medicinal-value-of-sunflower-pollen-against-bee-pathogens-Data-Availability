####analysis of Exp 9_ USA_WI_SUNFOWER POLLEN EXPERIMENT ####### 
###individual bees.... JONATHAN GIACOMINI.... 4 POLLEN DIETS EFFECTS ON CRITHIDIA COUNTS#
#1 - SUN - WI Sunflower Pollen
#2 - CSUN - China Sunflower Pollen
#3 - MIX- WI Mix POllen
#4 - KMIX - Koppert Mix Pollen 

#edits 2017.07.01: reorder letters for cld

##Edited EPY 2017.06.03 : use glmmTMB to fit model and drop1() to test predictor significance

#Edited 2017.06.04 to remove "screwed up" WI Mix pollen
#Edit 2017.06.11 for bar chart

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

#Remove WI MIX
Data<-subset(Data, !Pollen =="MIX")

#Change factor levels...
levels(Data$Pollen) ## "CSUN" "KMIX" "MIX"  "SUN" 
Data$Pollen<-factor(Data$Pollen, levels = c("KMIX", "MIX", "CSUN", "SUN"))


library(glmmTMB)
Data1<-na.omit(Data) #The following models don't work unless I omit the NAs


Data1<-droplevels(Data1)

#Poisson:
M.p<- glmmTMB(Count ~ Pollen   + wingmm + (1|ColonyID), 
              family = "poisson",
              data=Data1) 
summary(M.p)
# AIC      BIC   logLik deviance df.resid 
# 3057.6   3069.6  -1523.8   3047.6       75 
##Deviance >> df.resid .....---> overdispersed response variable

M.nb<-glmmTMB(Count ~ Pollen   + wingmm + (1|ColonyID), 
              family = "nbinom2",
              data=Data1) 
summary(M.nb)

AICctab(M.p, M.nb) 
# dAICc  df
# M.nb    0.0 6 
# M.p  2360.1 5 

#Stick with negative binomial then
drop1(M.nb, test = "Chisq")
# Df    AIC     LRT  Pr(>Chi)    
# <none>    697.17                      
# Pollen  2 715.71 22.5377 1.276e-05 ***
#   wingmm  1 695.28  0.1122    0.7376    
#wingmm not helpful; recode model without it, and
  #Bring back bees without wing measures
M.nb2<-update(M.nb, ~. -wingmm, data = Data)
drop1(M.nb2, test = "Chisq")
# Count ~ Pollen + (1 | ColonyID)
# Df    AIC    LRT  Pr(>Chi)    
# <none>    707.63                     
# Pollen  2 727.29 23.654 7.304e-06 ***

Lsm<-lsmeans(M.nb2, ~Pollen)
Pairs<-contrast(Lsm, "pairwise") 
Pairs
# contrast     estimate        SE df z.ratio p.value
# KMIX - CSUN  1.920254 0.3911556 NA   4.909  <.0001
# KMIX - SUN   1.665337 0.4219359 NA   3.947  0.0002
# CSUN - SUN  -0.254917 0.4239761 NA  -0.601  0.8194

##Sunflowers distinguish themselves

##Get letters to plot
cld(Lsm, Letters = letters, sort = FALSE)
#this will NOT rearrange the groups ... sort = TRUE will rearrange in ascending order

Lets<-cld(Lsm, Letters = letters, sort = FALSE)
##Edit 2017.06.27: try "sort" command
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
# Pollen   lsmean        SE df asymp.LCL asymp.UCL    Count      btlo     bthi .group
# 1   KMIX 4.502779 0.3967374 NA  3.725188  5.280370 90.26762 60.705936 134.2248      b
# 2   CSUN 2.582524 0.4074562 NA  1.783925  3.381124 13.23049  8.802785  19.8853      a
# 3    SUN 2.837441 0.4356649 NA  1.983554  3.691329 17.07203 11.042776  26.3932      a
Lsmeans.df<-Lsmeans.df.merged #overwrites old

#x-axis:
pollen.axis<-c("WF Mix","Sun (CN)", "Sun (USA)")
#y-axis:
#y-axis:
ylabel<- expression(italic(Crithidia)~count~
                      "(cells * 0.02"~mu~L^-1*")", sep="")  #check for encoding
Colors <- c("gray", "orange", "orange")

#ready to plot?
pd<-position_dodge(width=0.9)
p<-ggplot(Lsmeans.df, aes(x = Pollen, y = Count)) + 
  geom_bar(aes(fill = Pollen), color = "black", 
           stat = "identity", position = pd)+
  geom_errorbar(aes(ymin = btlo, ymax = bthi), 
                size = 1, width=0.4, color = "black", position = pd) +
  scale_x_discrete(labels=pollen.axis) +
  ylab(ylabel) + #y label
  xlab("Pollen")+
  theme(legend.position = "none")+
  theme(axis.line = element_line(size = 2),
        text = element_text (face = "bold"))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.2 * max(Lsmeans.df$bthi)))+
  scale_fill_manual(values = Colors)
p
#Add text
Panel.1d<-p + geom_text(aes(y = bthi + 0.1 * max(bthi), label = .group),
                        size = 8)
Panel.1d
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/NCSU_Crithidia/Expt9")
#ggsave("Exp9.china.usa.v4.bars.pdf", height = 5, width =5)



##################################################################################################
######### Binomial model for proportion of bees infected with Crithidia ###################
##################################################################################################
##################################################################################################


Data1$Infection<-Data1$Count
Data1$Infection[Data1$Infection>0] <- 1
Data1$Infection

aggregate(FUN = mean,Data1$Infection~Data1$Pollen)
# Data1$Pollen Data1$Infection
# 1         KMIX       0.9642857
# 2         CSUN       0.5925926
# 3          SUN       0.6800000

ProportionMod1<-glmmTMB(Infection~Pollen  + wingmm +
                          (1|ColonyID),  family = "binomial",
                        data=Data1) 

summary(ProportionMod1)


drop1(ProportionMod1, test = "Chisq")
# Single term deletions
# 
# Model:
#   Infection ~ Pollen + wingmm + (1 | ColonyID)
# Df    AIC     LRT  Pr(>Chi)    
# <none>    80.838                      
# Pollen  2 93.128 16.2898 0.0002902 ***
#   wingmm  1 79.126  0.2886 0.5911451     

ProportionMod2<-glmmTMB(Infection~Pollen + 
                          (1|ColonyID),  family = "binomial",
                        data=Data1) 

summary(ProportionMod2)
drop1(ProportionMod2, test = "Chisq")
# Single term deletions
# 
# Model:
#   Infection ~ Pollen + (1 | ColonyID)
# Df    AIC    LRT  Pr(>Chi)    
# <none>    79.126                     
# Pollen  2 91.201 16.075 0.0003231 ***
# 

#Source Ben Bolker function
source("glmmTMB.to.lsmeans.R")
library(lsmeans)
lsm.Prop<-lsmeans(ProportionMod2, ~Pollen)
Pairs<- contrast(lsm.Prop, "pairwise")
Pairs
# contrast      estimate        SE df z.ratio p.value
# KMIX - CSUN  3.4441767 1.1819334 NA   2.914  0.0100
# KMIX - SUN   2.9093591 1.1734857 NA   2.479  0.0351
# CSUN - SUN  -0.5348176 0.6456801 NA  -0.828  0.6854

LetsProp<-cld(lsm.Prop, Letters = letters)
LetsProp

# Pollen    lsmean        SE df  asymp.LCL asymp.UCL .group
# CSUN   0.2859578 0.7153674 NA -1.1161365  1.688052  a    
# SUN    0.8207754 0.7433168 NA -0.6360988  2.277650  a    
# KMIX   3.7301345 1.2258591 NA  1.3274949  6.132774   b  


#Use lsmeans to get mean counts
My.LsmeansProp <- lsmeans(ProportionMod2, ~Pollen) 
My.LsmeansProp
library(boot)
#inv.logit(x)
LsmeansProp.df<-as.data.frame(summary(lsmeans(ProportionMod2, ~Pollen)))
LsmeansProp.df$btmean<-inv.logit(LsmeansProp.df$lsmean)
LsmeansProp.df$btlo<-inv.logit(LsmeansProp.df$lsmean - LsmeansProp.df$SE) #lower bound
LsmeansProp.df$bthi<-inv.logit(LsmeansProp.df$lsmean + LsmeansProp.df$SE) #upper bound
View(LsmeansProp.df)
LsmeansProp.df

LetsProp2<- dplyr::select(LetsProp, Pollen, .group)
LetsProp2
library(plyr)
LsmeansProp.df.merged<- join(LsmeansProp.df, LetsProp2, by = "Pollen")
LsmeansProp.df.merged$.group<-gsub(pattern = " " , replacement="", x= LsmeansProp.df.merged$.group)
LsmeansProp.df.merged
# Pollen    lsmean        SE df  asymp.LCL asymp.UCL    btmean      btlo      bthi .group
# 1   KMIX 3.7301345 1.2258591 NA  1.3274949  6.132774 0.9765724 0.9244410 0.9930081      b
# 2   CSUN 0.2859578 0.7153674 NA -1.1161365  1.688052 0.5710063 0.3942673 0.7313191      a
# 3    SUN 0.8207754 0.7433168 NA -0.6360988  2.277650 0.6944009 0.5193550 0.8269398      a

#### PLOT- Proportion infected model ####
#x-axis:
pollen.axis<-c("WF Mix","SUN (PRC)","SUN (USA)")
#y-axis:
ylabel<- expression(italic(Crithidia)~prevalence, sep="")  #check for encoding
Colors <- c("gray", "orange", "orange")

#ready to plot?
pd<-position_dodge(width=0.8)
pX<-ggplot(LsmeansProp.df.merged, aes(x = Pollen, y = btmean)) +
  geom_bar(aes(fill = Pollen), color = "black", 
           stat = "identity", position = pd)+
  geom_errorbar(aes(ymin = btlo, ymax = bthi), 
                size = 1.25, width=0.6, color = "black", position = pd) +
  scale_x_discrete(labels=pollen.axis) +
  ylab(ylabel) + #y label
  xlab("Pollen")+
  theme(legend.position = "none")+
  theme(axis.line = element_line(size = 2),
        text = element_text (size = 20, face = "bold"))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1))+
  scale_fill_manual(values = Colors)
pX
PX2<-pX + geom_text(aes(y = bthi + 0.05 * max(bthi), label = .group),
                    size = 10)
PX2
PX3 <- PX2 + theme(text = element_text(size = 30, face="bold"))
PX3
PX4 <- PX3 + theme(axis.text.x = element_text(size=25), axis.text.y =element_text(size=25))
PX4

ggsave("Crithidia Prevalence.pdf", height = 9, width =9)

################# Composite figures #####################
#Use  plot_grid

#Must run other scripts or load workspaces:
#setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Analysis.Individual.Bee.and.Farm/Analyses/Individual Bees/Crithidia Load")
#"negative_binomial_individuals_v2.TMB.R"
#setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Analysis.Individual.Bee.and.Farm/Analyses/Individual Bees/Crithidia Load")
#negative_binomial_individuals_v3.barplot
#setwd("C:/Users/Evan/Google Drive/Project_Sunflower/NCSU_Crithidia/Exp2_JCRoulstonStrain_BuckSunKMix_2Week")
#"Exp2.buck.kopp.sun.TMB.v1a.R"


Panel.A2d<-plot_grid(Panel.1a, Ph2, Panel.1c, Panel.1d)
Panel.A2d

##Notes: Panels are 
#(A) Monofloral vs Mixed pollen
#(B) Individual survival
#(C) Effect of diet post-infection
#(D) Consistency of sunflower effect 
#(Supp) Reproducibility of effect... former panel C

#To-do:
##Add Panel Titles
##Remove x axis title
##Remove redundant y-axis titles

F1a.titled <- 
  Panel.1a + 
  ggtitle("Monofloral vs mixed pollen")+
  theme(axis.title.x = element_blank())
F1a.titled
F1b.titled<-Ph2 + ggtitle("Individual mortality")
F1b.titled
F1c.titled <- 
  Panel.1c + 
  ggtitle("Effect of diet post-infection")+
  theme(axis.title.x = element_blank()) #,
        #axis.title.y = element_blank())
F1c.titled
F1d.titled <- 
  Panel.1d + 
  ggtitle("Consistency of sunflower effect")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
F1d.titled
F1.all<- plot_grid(F1a.titled, F1b.titled, F1c.titled, F1d.titled,
                   nrow = 2, 
                   rel_widths = c(1, 1), 
                   align = "h",
                   labels = LETTERS, label_size = 20,
                   hjust = -0.4, vjust = 1.2)
F1.all

#where to put this??
# setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Manuscript Drafts/Figures")
# ggsave(F1.all, file = "Fig.1.4panel.v3.bars.hazards.pdf", height = 7, width = 7)



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


