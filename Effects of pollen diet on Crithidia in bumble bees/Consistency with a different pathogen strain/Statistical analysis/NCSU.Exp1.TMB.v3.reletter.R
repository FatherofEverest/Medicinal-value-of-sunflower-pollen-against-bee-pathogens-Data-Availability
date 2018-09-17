####analysis of POLLEN EXPERIMENT_1_JC Roulston Crithidia Strain ####### 
###individual bees.... JONATHAN GIACOMINI.... 2 POLLEN (Sun + Buck ) DIETS EFFECTS ON CRITHIDIA COUNTS#

###Edited by EPY 2017.06.03 : Models with glmmTMB
#Edit EPY 2017 06 11: barchart
#Edit 2017.07.01: chance lettering for cld

#Libraries:
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


###### Fetch data #######
#JG
data1<-read.csv("/Users/rg9403/Dropbox/Sunflower Pollen Experiments/Exp1/Sunflower Pollen Exp 1_Analysis 1 data.csv",header=TRUE)
#EPY
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/NCSU_Crithidia/Exp1_JCRoulstonStrain_BuckSun")
data1<-read.csv("Sunflower Pollen Exp 1_Analysis 1 data.csv",header=TRUE)

str(data1) 
# 'data.frame':	148 obs. of  6 variables:
#   $ BeeID    : int  1 2 3 4 5 6 7 8 9 10 ...
# $ ColonyID : Factor w/ 3 levels "JGA","JGB","JGC": 1 1 1 1 1 1 1 1 1 1 ...
# $ Date     : int  15336 15336 15336 15336 15336 15336 15336 15336 15336 15336 ...
# $ Treatment: Factor w/ 2 levels "BUCK","SUN": 2 1 2 1 2 1 2 1 2 1 ...
# $ Count    : num  0 0 0 78 0 0 0 1 0 0 ...
# $ WingSize : num  2.3 NA 2.3 2.1 1.95 2.1 2.05 2.1 2.4 2.2 ...

#Boxplots
car::Boxplot(Count~Treatment, data = data1)
 ##Low counts
d1.sum <- ddply(data1, ~Treatment, summarize, Median = median(Count))
d1.sum
# Treatment Median
# 1      BUCK      0
# 2       SUN      0

library(glmmTMB)
############## Negative Binomial Fit ##################
data1$Date<-factor(data1$Date) # Change Date to a factor
Data1<-rename(data1, Pollen = Treatment)

Fit1<-glmmTMB(Count ~ Pollen  + WingSize +(1|ColonyID)+(1|Date), 
               family = "nbinom2",
               data=Data1)
summary(Fit1)

# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  12.7445     3.5291   3.611 0.000305 ***
#   PollenSUN    -4.9868     0.9078  -5.493 3.95e-08 ***
#   WingSize     -4.7191     1.5634  -3.019 0.002540 ** 

drop1(Fit1)
##oops...
# Error in drop1.default(Fit1) : 
#   number of rows in use has changed: remove missing values?

#Refit without NA's
Fit.1a<-update(Fit1, data = na.omit(Data1))
drop1(Fit.1a, test = "Chisq")
# Count ~ Pollen + WingSize + (1 | ColonyID) + (1 | Date)
# Df    AIC    LRT  Pr(>Chi)    
# <none>      459.91                     
# Pollen    1 488.63 30.721 2.978e-08 ***
#   WingSize  1 468.08 10.171  0.001427 ** 

# We can knockout Date and colony ID and then compare models
Fit2<-update(Fit.1a, ~. - (1|Date))

Fit3<-update(Fit.1a, ~. - (1|ColonyID))

AICtab(Fit.1a, Fit2, Fit3)
# dAIC df
# Fit3   0.0  5 
# Fit2   1.6  5 
# Fit.1a 2.0  6 
 #Differences <2 AIC units. Let's keep Both Colony ID and Date as random effects
   #because they are part of the design

################### Plots #############
Lsm<-lsmeans(Fit.1a, ~Pollen)
Pairs<-contrast(Lsm, "pairwise") 
Pairs
# contrast   estimate        SE df z.ratio p.value
# BUCK - SUN 4.986755 0.9077968 NA   5.493  <.0001

##Sunflowers distinguish themselves

##Get letters to plot
Lets<-cld(Lsm, Letters = letters, sort = FALSE)
Lets
Lets$.group<-gsub(" ", "", Lets$.group) ##remove spaces

################## PLOT ON SCALE OF RESPONSE ##########
#Let's exponentiate original mean+SE and mean-SE
Lsmeans.df<-as.data.frame(summary(lsmeans(Fit.1a, ~Pollen)))
Lsmeans.df$Count<-exp(Lsmeans.df$lsmean)
Lsmeans.df$btlo<-exp(Lsmeans.df$lsmean - Lsmeans.df$SE) #lower bound
Lsmeans.df$bthi<-exp(Lsmeans.df$lsmean + Lsmeans.df$SE) #upper bound

Lsmeans.df

#Merge with letters
Lets2<- dplyr::select(Lets, Pollen, .group)
Lets2

Lsmeans.df.merged<- join(Lsmeans.df, Lets2, by = "Pollen")

Lsmeans.df.merged
# Pollen    lsmean       SE df asymp.LCL asymp.UCL      Count        btlo       bthi .group
# 1   BUCK  2.916269 0.586824 NA  1.766115  4.066423 18.4722401 10.27224007 33.2180374      b
# 2    SUN -2.070486 0.772864 NA -3.585271 -0.555700  0.1261245  0.05823029  0.2731808      a
Lsmeans.df<-Lsmeans.df.merged #overwrites old

#x-axis:
pollen.axis<-c("Buck", "Sun")
#y-axis:
ylabel<- expression(italic(Crithidia)~count~
 "(cells * 0.02"~mu~L^-1*")", sep="")  #check for encoding

#ready to plot?
pd<-position_dodge(width=0.9)
p<-ggplot(Lsmeans.df, aes(x = Pollen, y = Count)) + 
  #Option to sprinkle raw data... 
  geom_bar(aes(fill = Pollen), color = "black", 
           stat = "identity", position = pd)+
  geom_errorbar(aes(ymin = btlo, ymax = bthi), 
                size = 1, width=0.4, color = "black", position = pd) +
  scale_x_discrete(labels=c("Buckwheat", "Sunflower")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.2 * max(Lsmeans.df$bthi)))+
  scale_fill_manual(values = c("blue", "orange"))+
  ylab(ylabel) + #y label
  xlab("Pollen")+
  theme(legend.position = "none")+
  theme(axis.line = element_line(size = 2),
        text = element_text (face = "bold"))
p
#Add text
Panel.1d <- p + geom_text(aes(y = bthi + 0.1 * max(bthi), label = .group),
                          size = 8)
Panel.1d
# setwd("C:/Users/Evan/Google Drive/Project_Sunflower/NCSU_Crithidia/Exp1_JCRoulstonStrain_BuckSun")
# ggsave("NCSU.Exp1.v3.bars.pdf", height = 5, width =5)
# setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Manuscript Drafts/Figures")
# ggsave("NCSU.Exp1.v3.bars.pdf", height = 5, width =5)


