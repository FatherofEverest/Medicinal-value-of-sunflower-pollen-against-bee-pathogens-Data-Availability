####analysis of POLLEN EXPERIMENT_2_JC Roulston Crithidia Strain ####### 
###individual bees.... JONATHAN GIACOMINI.... 3 POLLEN (Sun + Buck + KMix) DIETS EFFECTS ON CRITHIDIA COUNTS#
### All bees were inoculated with Crithidia bombi on day 1 and day 2 (double dose)
### Bees were fed Kmix for first week and either switched to Buck or Sun, or remained on Kmix, for an additional week

#2017.07.01: reorder letters for cld
###EPY edits 2017.06.03: changed model to TMB
  ##notes: This uses only "Round2", there was another Round (Round1) that was dropped (why? Death?)
    #Round1 included a second inoculation @ 7 d p.i. ... also some deaths due to dehydration (ice storm prevented daily feed)

#EPY
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/NCSU_Crithidia/Exp2_JCRoulstonStrain_BuckSunKMix_2Week")

#Packages
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


####### Fetch data ######################
#D<-read.csv("/Users/rg9403/Dropbox/Sunflower Pollen Experiments/Sunflower Pollen Exp 2_Round2.csv",header=TRUE)
#str(D) 

#EPY:
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/NCSU_Crithidia/Exp2_JCRoulstonStrain_BuckSunKMix_2Week")
D<-read.csv("Sunflower Pollen Exp 2_Round2.csv")



str(D)
Boxplot(D$GutCount~D$Treatment) #cool

#Rename columns
D<-dplyr::rename(D, Pollen = Treatment)
D<-dplyr::rename(D, Count=GutCount)

############## Negative Binomial Fit ##################
Mod1<-glmmTMB(Count ~ Pollen  + WingSize +(1|ColonyID), 
                  family = "nbinom2",
                  data=D)
summary(Mod1)
drop1(Mod1, test = "Chisq")
# Count ~ Pollen + WingSize + (1 | ColonyID)
# Df    AIC     LRT  Pr(>Chi)    
# <none>      470.99                      
# Pollen    2 484.16 17.1772 0.0001862 ***
#   WingSize  1 470.60  1.6164 0.2035952    

### WingSize not sig (P = 0.2), drop it
Mod2<-update(Mod1, ~. - WingSize)
summary(Mod2)

drop1(Mod2, test = "Chisq")
# Df    AIC    LRT  Pr(>Chi)    
# <none>    470.60                     
# Pollen  2 482.25 15.643 0.0004009 ***


AICtab(Mod1, Mod2) #Not much difference 

#post hoc test 
summary(glht(Mod2,linfct=mcp(Pollen="Tukey")))
#Unhappy
# Error in modelparm.default(model, ...) : 
#   dimensions of coefficients and covariance matrix don't match
  #....hazards of work with brand new package

Lsm.mod<-lsmeans(Mod2, pairwise~Pollen)
Lsm.mod$contrasts
# contrast    estimate        SE df z.ratio p.value
# BUCK - KMIX 0.363062 0.5323870 NA   0.682  0.7740
# BUCK - SUN  2.456019 0.5290652 NA   4.642  <.0001
# KMIX - SUN  2.092957 0.5292723 NA   3.954  0.0002
# P value adjustment: tukey method for comparing a family of 3 estimates 
Pairs<-Lsm.mod$contrasts


###################################### PLOTS ####################################################

#Use lsmeans to get mean counts
My.Lsmeans <- lsmeans(Mod2, ~Pollen) 
My.Lsmeans
Lets<-cld(My.Lsmeans, Letters = letters, sort = FALSE)
Lets$.group<-gsub(" ", "", Lets$.group)

################## PLOT ON SCALE OF RESPONSE ##########
#Let's exponentiate original mean+SE and mean-SE
Lsmeans.df<-as.data.frame(summary(lsmeans(Mod2, ~Pollen)))
Lsmeans.df$Count<-exp(Lsmeans.df$lsmean)
Lsmeans.df$btlo<-exp(Lsmeans.df$lsmean - Lsmeans.df$SE) #lower bound
Lsmeans.df$bthi<-exp(Lsmeans.df$lsmean + Lsmeans.df$SE) #upper bound

Lsmeans.df

#Merge with letters
Lets2<- dplyr::select(Lets, Pollen, .group)
Lets2

Lsmeans.df.merged<- join(Lsmeans.df, Lets2, by = "Pollen")

Lsmeans.df.merged
# Pollen   lsmean        SE df asymp.LCL asymp.UCL     Count      btlo      bthi .group
# 1   BUCK 4.826813 0.3763089 NA  4.089261  5.564364 124.81251 85.670094 181.83897      B
# 2   KMIX 4.463751 0.3766000 NA  3.725628  5.201873  86.81251 59.569917 126.51371      B
# 3    SUN 2.370793 0.3718893 NA  1.641904  3.099683  10.70588  7.380962  15.52859     A 


Lsmeans.df<-Lsmeans.df.merged #overwrites old

#x-axis:
pollen.axis<-c( "Buck","WF Mix","Sun")
#y-axis:
ylabel<- expression(italic(Crithidia)~count~
                      "(cells * 0.02"~mu~L^-1*")", sep="")  #check for encoding
Colors <- c("blue", "gray", "orange")

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
Panel.1c<-p + geom_text(aes(y = bthi + 0.1 * max(bthi), label = .group),
                        size = 8)
Panel.1c
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/NCSU_Crithidia/Exp2_JCRoulstonStrain_BuckSunKMix_2Week")
# ggsave("Exp2_v3.barchart.pdf", height = 5, width =5)

