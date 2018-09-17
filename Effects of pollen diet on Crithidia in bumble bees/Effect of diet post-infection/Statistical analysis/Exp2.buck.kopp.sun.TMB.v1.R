####analysis of POLLEN EXPERIMENT_2_JC Roulston Crithidia Strain ####### 
###individual bees.... JONATHAN GIACOMINI.... 3 POLLEN (Sun + Buck + KMix) DIETS EFFECTS ON CRITHIDIA COUNTS#
### All bees were inoculated with Crithidia bombi on day 1 and day 2 (double dose)
### Bees were fed Kmix for first week and either switched to Buck or Sun, or remained on Kmix, for an additional week

data<-read.csv("/Users/rg9403/Dropbox/Sunflower Pollen Experiments/Sunflower Pollen Exp 2_Round2.csv",header=TRUE)
str(data) 

#EPY:
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/NCSU_Crithidia/Exp2_JCRoulstonStrain_BuckSunKMix_2Week")
D<-read.csv("Sunflower Pollen Exp 2_Round2.csv")

library(glmmTMB)
library(car)
library(bbmle) #Aic tables
library(lsmeans) #post-hoc tests with finite size correction
library(multcomp) #post-hoc tests
library(plyr)
library(dplyr)
source("glmmTMB.to.lsmeans.R") #helper function for TMB to talk to lsmeans

str(D)
Boxplot(D$GutCount~D$Treatment) #cool
D<-rename(D, Pollen = Treatment)
D<-rename(D, Count=GutCount)

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
  #hazards of work with brand new package

Lsm.mod<-lsmeans(Mod2, pairwise~Pollen)
Lsm.mod$contrasts
# contrast    estimate        SE df z.ratio p.value
# BUCK - KMIX 0.363062 0.5323870 NA   0.682  0.7740
# BUCK - SUN  2.456019 0.5290652 NA   4.642  <.0001
# KMIX - SUN  2.092957 0.5292723 NA   3.954  0.0002
# P value adjustment: tukey method for comparing a family of 3 estimates 
Pairs<-Lsm.mod$contrasts
#Letters for plot
Lets<-lsmeans::cld(lsmeans(Mod2, ~Pollen, Letters = c("a", "b", "c")))
Lets
Lets$.group<-gsub(" ", "", Lets$.group)
Lets  
#Why won't it use letters?

###################################### PLOTS ####################################################
library(lsmeans)
library(ggplot2)
#Use lsmeans to get mean counts
My.Lsmeans <- lsmeans(nbinom.TMB, ~Pollen) 
My.Lsmeans
cld(My.Lsmeans, Letters = LETTERS)


################## PLOT ON SCALE OF RESPONSE ##########
#Let's exponentiate original mean+SE and mean-SE
Lsmeans.df<-as.data.frame(summary(lsmeans(nbinom.TMB, ~Pollen)))
Lsmeans.df$Count<-exp(Lsmeans.df$lsmean)
Lsmeans.df$btlo<-exp(Lsmeans.df$lsmean - Lsmeans.df$SE) #lower bound
Lsmeans.df$bthi<-exp(Lsmeans.df$lsmean + Lsmeans.df$SE) #upper bound

Lsmeans.df
# Pollen    lsmean        SE df  asymp.LCL asymp.UCL      Count      btlo       bthi
# 1   Buck 4.6955379 0.2819003 NA 4.14302337  5.248052 109.457665 82.569266 145.102179
# 2    Mix 2.5084410 0.3030495 NA 1.91447487  3.102407  12.285762  9.073803  16.634695
# 3   Rape 3.7538852 0.2889668 NA 3.18752063  4.320250  42.686604 31.973849  56.988641
# 4    Sun 0.7166937 0.3197605 NA 0.08997472  1.343413   2.047652  1.487257   2.819203

#Merge with letters
Lets2<- dplyr::select(Lets, Pollen, .group)
Lets2

Lsmeans.df.merged<- join(Lsmeans.df, Lets2, by = "Pollen")

Lsmeans.df.merged
Lsmeans.df.merged$.group<-gsub(pattern = " " , replacement="", x= Lsmeans.df.merged$.group)
Lsmeans.df.merged

Lsmeans.df<-Lsmeans.df.merged
#y-axis:

#x-axis:
pollen.axis<-c( "Buck","Rape","Mix", "Sun")
#y-axis:
ylabel<- expression(paste("Parasite load (ln (cells * 0.02 ", mu, L^-1, "))", sep="") ) #check for encoding

#ready to plot?
p<-ggplot(Lsmeans.df, aes(x = Pollen, y = Count)) + 
  #Option to sprinkle raw data... 
  geom_point(data = Data1, aes( x = Pollen, y = Count), position = "jitter", 
             color = "blue", shape = 25, size = 0.5, alpha = 0.4)+
  geom_pointrange(aes(ymin = btlo, ymax = bthi), size = 1, color = "black") +
  scale_x_discrete(labels=pollen.axis) +
  ylab(ylabel) + #y label
  xlab("Pollen")
p
#Add text
p + geom_text(aes(y = bthi + 0.2 * max(bthi), label = .group),
              size = 8)
ggsave("counts_individ_v2.TMB.pdf", height = 5, width =5)





######## PLOT ########
library(lsmeans)
My.Lsmeans <- lsmeans(Mod2, ~Treatment) 
My.Lsmeans
# Treatment   lsmean  SE       df asymp.LCL asymp.UCL
# BUCK      4.826868 0.3763100 NA  4.089314  5.564422
# KMIX      4.463766 0.6494017 NA  3.190962  5.736570
# SUN       2.370775 0.6538854 NA  1.089183  3.652367

library(ggplot2)
Lsmeans.df<-summary(My.Lsmeans)
str(Lsmeans.df)
View(Lsmeans.df)
#label for x axis
pollen.axis<-c( "Buckwheat","K-Mix", "Sunflower")

plot0<- ggplot(Lsmeans.df, aes(x=Treatment, y=lsmean)) + 
  geom_bar(position=position_dodge(), stat="identity", colour="black" ) +
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE),
                size=1, width=0.2,  # Width of the error bars
                position=position_dodge(0.9))
plot0 # Looks ok, let's refine it now
plot01<-plot0 + theme_bw() #remove colored background
ylabel<- expression(bold(paste("Parasite load (ln (cells * 0.02 ", L^-1, "))", sep="") )) #check for encoding
plot02<-plot01 + ylab(ylabel) #y label
plot02
#relabel x axis:
plot05<-plot02 + xlab("Pollen type") + # x axis label
  scale_x_discrete(breaks=levels(Lsmeans.df$Treatment) , labels=pollen.axis)
plot05
#Larger font:
plot06<-plot05 + theme(text=element_text (size=25, face="bold") )
plot06
#Remove gridlines
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
fig1<-Justified
fig1

