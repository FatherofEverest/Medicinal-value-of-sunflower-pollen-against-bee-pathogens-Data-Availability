####2016.03#########
#analysis of sunflower microcolony performance

#####Egg production################
####Source data: use summary created by LS Adler
#EPY added back colonies from which new adults emerged -- 
#colonies 3, 11, 23, 25, 33
#coded as "1" for pupal production

rm(list=ls()) #clear memory
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/jess stats/colony.performance")

#Libraries:
##Load packages 
library(MASS)
library(glmmTMB)
library(lme4)
library(car)
library(bbmle) #Aic tables
library(lsmeans) #post-hoc tests with finite size correction
library(multcomp) #post-hoc tests
library(plyr)
library(dplyr)
source("glmmTMB.to.lsmeans.R") #helper function for TMB to talk to lsmeans
library(ggplot2)
library(cowplot) #multi panel plots


Perf0<- read.csv('C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/jess stats/colony.performance/performance_sum_SAS_CALLOWS.ADDED.csv',
                 na.strings=".")
Bees<-read.csv("C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/Pollen consumption trials/updated survival.csv")
#use "bees" to get a column for dimorphism
Data<-Bees
Data$wingmm<-(13/20)*Data$marginal.cell #converts wing measurement from ocular units to mm
#first make data frame with max and min size for each bee in the colony
#example:
Data$Microcolony<-Data$MC #renaming
library(plyr)
DM<-ddply(Data,~Microcolony, summarise,bigbee=max(wingmm),smallbee=min(wingmm))
str(DM) #good, data frame
#calculate colony dimorphism: ratio of largest to smallest bee minus 1
DM$dimorph<-(DM$bigbee/DM$smallbee) -1
DM
#merge this back to original data:
Perf0$Microcolony<-Perf0$microcolony
Data.DM<-join(x=Perf0, y=DM, by=c("Microcolony"), type="left")
Data.DM$Microcolony
Perf0<-Data.DM

#familiarize with the Performance data
str(Perf0)
#exclude microcolonies 12, 16 (uninfected treatment, found infected), 
#also exclude 50, 54 (died within first week)
Perf02<-subset(Perf0, !Microcolony=="50" & !Microcolony=="54" & !Microcolony=="12" & !Microcolony=="16")
Perf02$Microcolony<-as.factor(Perf02$Microcolony) #changes from numeric to factor
str(Perf02$Microcolony) #check that it worked
str(Perf02)
str(Perf02$sourcecol) #good , 4 levels
library(dplyr)
Perf02$sourcecol
#View(Perf02)
#rename(Perf02, Colony=sourcecol) #not working
Perf02$Colony<-Perf02$sourcecol
str(Perf02) #duplicated, renamed the column
Perf03<-Perf02
str(Perf03$infect) #good , 2 levels
Perf03$Infection<-Perf03$infect
Perf03$Infection<-relevel(Perf03$Infection, ref="U")
#reference level now "U"
str(Perf03) #renamed the column

Perf03$Pollen<-Perf03$pollen

#check eggnum
Perf03$eggnum
hist(Perf03$eggnum)
Buck<-subset(Perf03, Pollen=="B")
Sun<-subset(Perf03, Pollen=="S")
hist (Buck$eggnum)
hist(Sun$eggnum)
mean(Buck$eggnum) #<2 eggs
mean(Sun$eggnum) #2.5 eggs


library(car)
library(plyr)
?summarize #to summarize by a grouping variable
avg.eggs<-ddply(Perf03,~Pollen * Infection, summarise,Production=mean(eggnum))
avg.eggs
ggplot(avg.eggs, aes(x=Infection, y=Production))+
  geom_bar(stat="identity", color="black") + facet_grid(Pollen~.)
#certainly looks like a big interaction there!
Boxplot(eggnum~Pollen*Infection, data=Perf03)
#interesting-- infection INCREASED egg production in sunflower
#but CRATERED egg production in buckwheat!!

###MODELING STARTS HERE#####

#ready to run Poisson model with Colony as random effect?
 
#Set contrasts:
contrasts = c("contr.sum", "contr.poly")

egg3<-glmer( eggnum~ Infection * Pollen  + (1|Colony), data=Perf03,
             nAGQ = 25,
             family=poisson)
summary(egg3)
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)         0.85552    0.30960   2.763  0.00572 ** 
#   InfectionI         -1.12580    0.28059  -4.012 6.01e-05 ***
#   PollenS            -0.07778    0.19819  -0.392  0.69471    
# InfectionI:PollenS  1.10600    0.34401   3.215  0.00130 ** 
Anova(egg3, Type = 3)
# Chisq Df Pr(>Chisq)   
# Infection         5.7713  1   0.016290 * 
#   Pollen            3.1946  1   0.073883 . 
# Infection:Pollen 10.3366  1   0.001304 **

drop1(egg3, test = "Chisq")
# Df    AIC    LRT  Pr(Chi)    
# <none>              240.09                    
# Infection:Pollen  1 249.11 11.017 0.000903 ***

library(lsmeans)
margmeans<-summary(lsmeans(object=egg3, ~Infection*Pollen))
margmeans2<-summary(lsmeans(object=egg3, ~Infection*Pollen), type="response")
margmeans2
# Infection Pollen      rate        SE df asymp.LCL asymp.UCL
# U         B      2.3525919 0.7283665 NA 1.2823694  4.315986
# I         B      0.7631606 0.2809681 NA 0.3708808  1.570354
# U         S      2.1765380 0.6724155 NA 1.1879486  3.987814
# I         S      2.1338562 0.6605826 NA 1.1632064  3.914475

margmeans
margmeans$btmean<-exp(margmeans$lsmean)
margmeans$rawhi<-margmeans$lsmean + margmeans$SE
margmeans$rawlo<-margmeans$lsmean - margmeans$SE
margmeans$bthi<-exp(margmeans$rawhi)
margmeans$btlo<-exp(margmeans$rawlo)
margmeans$btmean #matches margmeans2, good
margmeans
#I wanted to check that I could
#get the back-transformed results by hand by exponentiating


#########plotting
library(lsmeans)
margmeans2
str(margmeans2)
margmeans2<-as.data.frame(margmeans2)
View(margmeans2)
My.Lsmeans <- margmeans2
My.Lsmeans

#######PLOTTING : egg number#############
#plotting marginal means:
library(ggplot2)
#make a data frame from the lsmeans object
Lsmeans.df<-My.Lsmeans
Lsmeans.df$lsmean<-Lsmeans.df$rate
str(Lsmeans.df)
View(Lsmeans.df) #???
levels(Lsmeans.df$Nectar.conc)

#below plot will give SYMMETRIC error bars:
#plot0<- ggplot(Lsmeans.df, aes(x=Pollen, y=lsmean, fill=Infection)) + 
 # geom_bar(position=position_dodge(), stat="identity", colour="black", size=1.5 ) + #color is for OUTLINES
 # geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), #adds error bars at 2x thickness
  #              size=2, width=0.2,  # Width of the error bars
  #              position=position_dodge(0.9))+ #useful command for line charts to offset error bars, vestigial here?
   #             theme_classic(base_size=25) #this theme has no background, default base_size is 12 point font


#Use this plot for manually exponentiated mean+SE from scale of linear predictor(asymmetric) 
plot0<- ggplot(margmeans, aes(x=Pollen, y=btmean, fill=Infection)) + 
  geom_bar(position=position_dodge(), stat="identity", colour="black", size=1.5 ) + #color is for OUTLINES
  geom_errorbar(aes(ymin=btlo, ymax=bthi), #adds error bars at 2x thickness
               size=2, width=0.2,  # Width of the error bars
               position=position_dodge(0.9))+ #useful command for line charts to offset error bars, vestigial here?
              theme_classic(base_size=25) #this theme has no background, default base_size is 12 point font
  

plot0 #Solid base plot
p01<-plot0 +   theme(legend.key = element_rect(fill = 'black', size=.1)) #black outline, thinner lines
p01 
p01big<-p01 + theme(legend.key.size=unit(1.5, "cm"))

plot02<-p01big + ylab("Number of eggs") + xlab("Pollen") #y label
plot03<-plot02 + theme(text=element_text(face="bold")) + theme(line=element_line(size=2))
Cbombi<-expression(paste(italic("Crithidia")))

p4<-plot03 + scale_fill_manual(name=Cbombi, #use object created above
                      values=c("white", "gray"), #colors for each series
                      labels=c("Uninfected", "Infected")) #Legend: full treatment names
#remove space below bars
p5<-p4 + coord_cartesian(ylim = c(0, 3.2))+
  scale_y_continuous(expand = c(0,0))
p5
#Pollen: full pollen name
#label for x axis
pollen.axis<-c("Buckwheat", "Sunflower")
p6<-p5 + scale_x_discrete(labels=pollen.axis)
p6
p6 + legend_theme
fig1<-p6

#give directory for output file
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/jess stats/colony.performance")

#export to PDF
#dev.new()
# pdf("eggfigure_manual_SE.pdf",height=7,width=7,paper='special') 
# fig1
# dev.off()

#for eps

#setEPS()
# postscript("eggfigure_manual_SE.eps", horizontal = FALSE, onefile = FALSE, 
#            paper = "special", height=7,width=7)
# fig1
# dev.off()


