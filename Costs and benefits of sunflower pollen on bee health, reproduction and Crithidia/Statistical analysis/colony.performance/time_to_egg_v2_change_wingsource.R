###TIME TO EGG LAYING####

#Data from J Leslie sunflower/buckwheat  crithidia microcolony experiment
#Microcolony	-- experimental unit of 5 adult bees
#Innoculationdate	-- date on which colony inoculated; corresponds to colony of origin
#Dateegg1 -- calendar date of egg laying	
#Time2egg	-- number of days after inoculation before first egg found
#Laideggs-- whether or not eggs were laid during the experiment	
#Eggbi	-- binary response , 0 = no egg, 1= laid egg
#Numbereggs--no data (number of eggs laid before dissection time 
#Infection: treatment with 2 levels, I=inoculated with Crithidia-containing gut extract, U=sham inoculated
#Pollen : pollen type that microcolonies were fed during experiment. B=buckwheat, S=sunflower


rm(list=ls())
#R packages to install, load
#install.packages("coxme") #do this for all packages that you need
library(coxme)
library(car)
library(plyr)

#working directory
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/jess stats/colony.performance")
EGG<-read.csv("updated preformance.csv")
Bees<-read.csv("C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/Pollen consumption trials/updated survival.csv")

summary(EGG)
#exclude microcolonies 50 and 54-- these died in first week post inoculation
Egg1<-subset(EGG, !Microcolony=="50" & !Microcolony=="54")
length(Egg1$Microcolony)
#exclude colonies 12 and 16, these were uninfected treatments that had infection.
Egg2<-subset(Egg1, !Microcolony=="12" & !Microcolony=="16")
length(Egg2$Microcolony)
#left with 76 microcolonies
plot(Time2egg~Pollen + Infection, data=Egg2) #lower for sunflower

#use data from individual bees to calculate size dimorphism of each microcolony
#bigbee	: wing size of largest bee in microcolony
#smallbee: : wing size of smallest bee in microcolony
#dimorph : measure of size dimorphism -- defined as  
# dimorph = (bigbee/smallbee) - 1
#calculate colony size dimorphism for use as covariate:
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
Data.DM<-join(x=Egg2, y=DM, by=c("Microcolony"), type="left")
Data.DM
length(Data.DM$Microcolony)#76
length(Data$Microcolony)#76, good
Data<-Data.DM

Measured<-subset(DM, !dimorph=="NA")
length(Measured$Microcolony) #lost 9 microcolonies-- too many!


#new column for "Treatment", combining infection treatment and pollen treatment
Data$Treatment <- paste (Data$Infection, Data$Pollen, sep=".")
Data$Treatment 
str(Data$Treatment)
levels(factor(Data$Treatment)) #good, 4 levels

Data<-transform(Data, Group=Treatment, Death_Days=Time2egg, Dead=Eggbi, Inoculation.date=Innoculationdate)

str(Data$Group)
####REORDER factor levels####
#U.B., U.S., I.B., I.S.
#now: IS, IB, UB, US
#change to 3, 4, 2, 1
Data$Group= factor(Data$Group,levels(Data$Group)[c(3,4,1,2)])
str(Data$Group)
summary(Data$Group)
str(Data$Group)

Data<-droplevels(Data) #to avoid errors caused by empty factor levels
Data$Microcolony
summary(Data)
summary(Data$Death_Days) #max of 42 days
#only 2 colonies did not lay eggs

Data$Microcolony<-factor(Data$Microcolony)
D1<-subset(Data, !dimorph=="NA")
Nested.Plus<-coxme(Surv(Death_Days, Dead) ~ 
                     Infection*Pollen +  dimorph + (1|Inoculation.date) , data=D1 )
library(car)
summary(Nested.Plus)
Anova(Nested.Plus)

#eliminate dimorphism
NP1<-update(Nested.Plus, ~. - dimorph)
summary(NP1)

#Now we can add back the colonies without dimorphism measure:
NP1a<-coxme(Surv(Death_Days, Dead) ~ 
              Infection*Pollen +  (1|Inoculation.date) , data=Data )
summary(NP1a)

#eliminate interaction
NP2<-update(NP1a, ~. - Infection:Pollen )
summary(NP2)



#Test the random effect:
#Remove random effect; advise against this because it is part of the design!!!
FixedFit<-coxph(Surv(Death_Days, Dead) ~ Infection + Pollen  , data=Data )
anova(NP2, FixedFit)
#Log-likelihood significantly higher (less negative) with the random effect; 
#KEEP the random effect

Final<-NP2
summary(NP2)
Anova(NP2)
#interpret coefficients (see links for info):
#http://stats.stackexchange.com/questions/6026/how-do-i-interpret-expb-in-cox-regression
#http://www.stat.nus.edu.sg/~stachenz/ST3242Notes4.pdf
#coefficients correspond to % increase in log hazard rate
#exp(coef) = exponentiated coefficients (ie e^(coef)) 
#reflect death hazard ratios relative to ref level (Uninfected, Buckwheat)
#at any particular point in time



#######################PLOTTING##########################
####BASE PLOT###########
Bee.Survival <- Surv(Data$Death_Days, Data$Dead)
#create the Kaplan Meier object:
Bee.Survival.Treatment <- survfit(Bee.Survival ~ Group, data=Data)
#Bee.Survival.Treatment <- survfit(Bee.Survival ~ Nec_Treat +Pol_Treat, data=Data)
summary(Bee.Survival.Treatment)
plot.new

#export plot
#choose source file loaction for wd
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/jess stats/colony.performance")
#uncomment for pdf:
#pdf("fig3_eggtime.pdf",height=7,width=7,paper='special') 
#setEPS() #uncomment this and next 2 lines for eps 
#postscript("fig2_survival.eps", horizontal = FALSE, onefile = FALSE, 
#paper = "special", height = 7, width = 7)

#Fuss with margins to remove excess white space:
par(mar=c(5, 5.5, 1, 1) + 0.2) #bottom, left, top, right
par(mgp=c(3,1,0))
plot(
  frame=FALSE,
  Bee.Survival.Treatment,fun="event",
  conf.int = FALSE,
  xlab=expression(bold("Time (d)")), 
  ylab=expression(bold("Proportion with egg")), 
  xlim = c(0,45), 
  ylim = c(0.0,1.02) ,
  col = c("black","black","black","black"),
  lty = c(2,2,1,1), #broken lines for uninfected, solid for infected
  pch=c(1,0,16,15), cex=2, #symbols to use
  mark.time=c(1,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53), # which timepoints to mark 
  mark=c(1,0,16,15), #again, which symbols
  lwd = 2,
  cex.axis = 1.75,
  cex.lab = 2.25
)
box(bty="l", lwd=2)
#surv.title<-expression(bold("Diet Treatment"))
#legend(x=2.0,y=0.75,  xjust=0.5, legend="", title = surv.title, cex=2, bty="n")
#this first "legend has just the title, in larger font (cex)
legend("bottomright",bty="n", #removes box around legend
       #x=1, y=0.7,       #for bottom left legend
       title = "", #no title 
       c("Uninfected, Buckweat", "Uninfected, Sunflower", "Infected, Buckweat", "Infected, Sunflower"),
       lwd = 2, cex = 1.5,
       col = c("black","black","black","black"),
       lty = c(2,2,1,1),
       pch=c(1,0,16,15), #symbols to use
)
#dev.off()





