####analysis of POLLEN EXPERIMENT ####### 
  ##First experiment, winter 2014-2015 at UMass Amherst ##
  ##Used Chinese (Huading Wax Industry) single and mixed pollens

##Updated 2017.06.03: Use glmmTMB instead of glmmADMB

###individual bees.... JONATHAN GIACOMINI.... 4 POLLEN DIETS EFFECTS ON CRITHIDIA COUNTS#
library(ggplot2)
library(cowplot)
library(glmmADMB)
library(glmmTMB)
library(lme4)
library(plyr)
library(dplyr)
library(lsmeans)


####Final model -> nbinom.TMB

#set working directory
rm(list=ls()) #clear memory

#EPY:
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Analysis.Individual.Bee.and.Farm/Analyses/Individual Bees/Crithidia Load")

#read in Jon's data from somewhere
#EPY:
data<-read.csv("Analysis_1_Crithidia load.csv",header=TRUE)
str(data)
data<-read.csv("/Users/rg9403/Google Drive/Project_Sunflower/Analysis.Individual.Bee.and.Farm/Analyses/Individual Bees/Crithidia Load/Analysis_1_Crithidia load.csv", header=TRUE)
str(data)

Data1<-data
Data2<-na.omit(Data1)
##Change levels of factor
str(Data1$Treatment)
Data1$Treatment <- factor(Data1$Treatment, levels = c("Buck", "Rape", "Mix", "Sun"))
levels(Data1$Treatment)

Data1<-rename(Data1, c("Treatment" = "Pollen"))
names(Data1)
Data2<-na.omit(Data1)
###### Refit with glmmTMB ###
nbinom.TMB<-glmmTMB(Count ~ Pollen  + WingSize +(1|ColonyID)  + (1|InocDate), 
                  family = "nbinom2",
                  data=Data2) 
#Use drop1() because car::Anova not supported
drop1(nbinom.TMB, test = "Chisq")
# Df    AIC     LRT Pr(>Chi)    
# <none>    1765.4                     
# Pollen  3 1870.6 111.239  < 2e-16 ***
#   wingmm  1 1770.0   6.644  0.00995 ** 


#Pairwise comparisons

#Source Ben Bolker function
source("glmmTMB.to.lsmeans.R")

lsm.TMB<-lsmeans(nbinom.TMB, ~Pollen)
Pairs<- contrast(lsm.TMB)
Pairs
# contrast      estimate        SE df z.ratio p.value
# Buck effect  1.7768984 0.1903050 NA   9.337  <.0001
# Mix effect  -0.4101984 0.2012969 NA  -2.038  0.0416
# Rape effect  0.8352457 0.1915766 NA   4.360  <.0001
# Sun effect  -2.2019457 0.2204260 NA  -9.989  <.0001
#P value adjustment: fdr method for 4 tests 

Lets<-cld(lsm.TMB, Letters = letters)
Lets
# Pollen    lsmean        SE df  asymp.LCL asymp.UCL .group
# Sun    0.7166937 0.3197605 NA 0.08997472  1.343413  a    
# Mix    2.5084410 0.3030495 NA 1.91447487  3.102407   b   
# Rape   3.7538852 0.2889668 NA 3.18752063  4.320250    c  
# Buck   4.6955379 0.2819003 NA 4.14302337  5.248052     d 

#compare to asking for response scale initially:
lsmeans(nbinom.TMB, ~Pollen, type = "response")
#doesn't work with TMB, it seems

###################################### PLOTS ####################################################
library(lsmeans)
library(ggplot2)
#Use lsmeans to get mean counts
My.Lsmeans <- lsmeans(nbinom.TMB, ~Pollen) 
My.Lsmeans

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

Lsmeans.df$Crithidiacells<-Lsmeans.df$Count * 50
Lsmeans.df$CrithidiacellsLowSe<-Lsmeans.df$btlo * 50
Lsmeans.df$CrithidiacellsHiSe<-Lsmeans.df$bthi * 50
Lsmeans.df
#Merge with letters
Lets2<- dplyr::select(Lets, Pollen, .group)
Lets2

Lsmeans.df.merged<- join(Lsmeans.df, Lets2, by = "Pollen")

Lsmeans.df.merged
Lsmeans.df.merged$.group<-gsub(pattern = " " , replacement="", x= Lsmeans.df.merged$.group)
Lsmeans.df.merged
# Pollen    lsmean        SE df  asymp.LCL asymp.UCL      Count      btlo       bthi Crithidiacells CrithidiacellsLowSe CrithidiacellsHiSe .group
# 1   Buck 4.6955379 0.2819003 NA 4.14302337  5.248052 109.457665 82.569266 145.102179      5472.8833          4128.46331          7255.1089      d
# 2   Rape 3.7538852 0.2889668 NA 3.18752063  4.320250  42.686604 31.973849  56.988641      2134.3302          1598.69246          2849.4320      c
# 3    Mix 2.5084410 0.3030495 NA 1.91447487  3.102407  12.285762  9.073803  16.634695       614.2881           453.69017           831.7347      b
# 4    Sun 0.7166937 0.3197605 NA 0.08997472  1.343413   2.047652  1.487257   2.819203       102.3826            74.36283           140.9601      a

Lsmeans.df<-Lsmeans.df.merged
#y-axis:

#x-axis:
pollen.axis<-c( "Buck","Rape","Mix", "Sun")
#y-axis:
ylabel<- expression(paste("Parasite load (cells * 0.02 ", mu, L^-1, "))", sep="") ) #check for encoding

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
Panel.1a<-p + geom_text(aes(y = bthi + 0.2 * max(bthi), label = .group),
              size = 8)
#ggsave("counts_individ_v2.TMB.pdf", height = 5, width =5)















