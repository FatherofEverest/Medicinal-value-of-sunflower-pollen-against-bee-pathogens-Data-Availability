#Analysis of Crithidia counts in microcolonies fed different pollens
#Jessica Leslie honors thesis 2015
#updated 2016.02.29 to remove "BeeID" random effect 
#(model was over-specified)
#script updated 2016.02.16 to use updated file "updated survival"
#mistakes where there were 2 rows for a single bee are now removed
#some wing measures have been added



#Adler notes from SAS file
#Variables:
#a microcolony expt in B impatiens with C bombi 
#conducted by Jess Leslie as an REU in the summer of 2015. 
#Each microcolony started with 5 workers from a colony of origin. 
# 2 infection treatments (U(ninfected) and I(nfected) 
#crossed  with two pollen diets: S(unflower, Helianthus annuus) or B(uckwheat, Fagopyrum - check species). 
#We had 5 microcolonies per treatment combination per colony of origin, 
#for a total of (right now) 20 colonies per colony of origin. 
#    BeeID is a roman number from I to V for each of the 5 workers in original microcolony. 
#Microcol is the microcolony number,
#infection and pollen are the treatments described above. 
# these are crithidia count data are in cells/uL (extrapolation from count of 0.02uL gut extract)
# Died is yes/no for whether it died during the trial; 
# Bees that died during trial were NOT DISSECTED for assessment of Crithidia
#note-- dead bees will be used for survival analysis). 
#Deathdate is the date of death. 
#crithidia is the Crithidia count, measured as 
#number of motile cells per 0.02 uL gut extract, multiplied by 50 for cells/uL
# Wingocun is the radial cell length of the right forewing in ocular units. For all the bees in this expt, 20 ou = 13 mm*/
# wingmm = (13/20)*wingocun

#variables in "countsforJMP"
#BeeID: roman numeral I-V,  arbitrarily assigned 	to each of the 5 bees in the microlony
#Microcol	: microcolony number-- each microcolony consisted of 5 adult bees
#colony: bee's natal colony; 20 microcolonies built from each colony
#Each colony's microcolonies were inoculated on the same day
#Died	: whether bee died during experiment; bees that died were nota assessed for Crithidia count
#Infection : Crithidia treatment-- two levels -- (U(ninfected) and I(nfected) 
#Pollen: two pollen diets: S(unflower, Helianthus annuus) or B(uckwheat, Fagopyrum - check species). 
#_TYPE_ : ???WHAT IS THIS
#_FREQ_ : ???WHAT IS THIS
#wingmm	: radial cell length of the right forewing in mm
#daysinfect : number of days between inoculatino and dissection, which occurred set # of days after first egg-laying
#totcrith : number of motile cells per 0.02 uL gut extract, multiplied by 50 for cells/uL
#bigbee	: wing size of largest bee in microcolony
#smallbee: : wing size of smallest bee in microcolony
#dimorph : measure of size dimorphism -- defined as  
# dimorph = (bigbee/smallbee) - 1

rm(list = ls())
#load packages
#install.packages("R2admb") #run if need to install, takes some time esp. next line
#install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos",getOption("repos")), type="source")
#run "install.packages("packagename") if you need to install-- remember the quotation marks!
library("glmmADMB") #for mixed model (glmmadmb)
library(lme4) #alternative mixed model, function "glmer"
library("ggplot2") #for plotting
library("car") #chi-square tests (Anova)
library("lsmeans")


#set working directory
rm(list=ls())
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/jess stats/Counts/R_analyses")
Data<-read.csv("C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/Pollen consumption trials/updated survival.csv")#rename some columns
View(Data)

#check dataset
str(Data)

#eliminate a priori microcolonies 12, 16 (uninfected treatment found infected)
# and 50, 54 (died within first week, did not collect data)
Data<- subset(Data, !MC==12 & !MC==16 & !MC==50 & !MC==54)
summary(Data$MC)

#rename some columns
Data$Microcolony<-Data$MC
Data$Microcol<-factor(Data$Microcol) #converts to factor
Data$Microcol #76 levels, good
Data$Colony<-Data$innoc..Date #inoculation date identifies which colony was used
Data$Pollen<-Data$pollen.treat
Data$daysinfect<-Data$days.infect
Data$totcrith<-Data$crithidia.count
#Data$totcrith<-Data$crithidia.count #converts count to cells/uL (option)
Data$wingmm<-(13/20)*Data$marginal.cell #converts wing measurement from ocular units to mm

Data
str(Data)

#now calculate size dimorphism for each colony
library(plyr)
#first make data frame with max and min size for each bee in the colony
#example:
#ddply(dt,~group,summarise,mean=mean(age),sd=sd(age))
DM<-ddply(Data,~Microcol, summarise,bigbee=max(wingmm),smallbee=min(wingmm))
str(DM) #good, data frame
#calculate colony dimorphism: ratio of largest to smallest bee minus 1
DM$dimorph<-(DM$bigbee/DM$smallbee) -1
DM
#merge this back to original data:
Data.DM<-join(x=Data, y=DM, by=c("Microcol"), type="left")
Data.DM
View(Data.DM)
length(Data$Microcol)#380
length(Data.DM$Microcol)#380, good
#use only infected bees in counts analysis:
Infected<-subset(Data.DM, !infect.treat=="U")
Infected$infect.treat #some problems here, "I" recognized as 2 levels
Infected<-droplevels(Infected)
str(Infected$infect.treat)
Infected$infect.treat
library(plyr)
Infected$infect.treat<-revalue(Infected$infect.treat, c("I "="I"))#merge the 2 levels
str(Infected$infect.treat) #good, down to 1 level
#eliminate bees without counts
Counted<-subset(Infected, !totcrith=="NA")


str(Counted$Pollen) #2 levels
Data<-Counted #renames
Data$BeeID<-Data$Bee.ID #renaming
Data$colony<-Data$Colony #renaming
Data<-droplevels(Data) #remove empty factor levels
############     MODELING        #######
#Start with: $totcrith as response variable
#$pollen as fixed effect, 
#Covariates of $daysinfect, $dimorph
#random effects: (blocking term accounting for non-independence of bees within microcolony within colony)
# + 1( colony|Microcol|BeeID )
#A model of nested random e???ects (block within site)
#would be 1|site/block; 
#a model of crossed random e???ects (block and year) would be:
#(1|block)+(1|year).

#Compare Poisson with zero inflation, regular negative binomial, and negative binomial with zero inflation
#Example:
#fit_zipoiss <- glmmadmb(NCalls ~ (FoodTreatment + ArrivalTime) * SexParent +
   #offset(log(BroodSize)) + (1 | Nest), data = Owls, zeroInflation = TRUE,
    #family = "poisson") 
#first try poisson model:
fitPoiss<- glmmadmb(totcrith ~ Pollen  + #fixed
                      daysinfect +dimorph + 1|Microcol, 
                    family = "poisson",
                    data=Data)  #throwing errors with poisson fit, no good
#try zero inflation with negative binomial:
zinfl<- glmmadmb(totcrith ~ Pollen + (1 | Microcol), data=Data,
                 family = "nbinom",zeroInflation = TRUE )
#this ran OK, but I don't think zero inflation model is appropriate
#because we experimentally inoculated all bees

#now negative binomial
#Start with full nesting
Winged<-subset(Data, !wingmm=="NA" ) #remove some NA's
Dimorphed<-subset(Winged, !dimorph=="NA") #remove some NA's

nbinom0<-glmmadmb(totcrith ~ Pollen  +  wingmm + dimorph + (1 |Colony/ Microcol), 
                  family = "nbinom",
                  data=Dimorphed)  #got an error, let's try to simplify
summary(nbinom0) #remove dimorphism
nbinom01<-glmmadmb(totcrith ~ Pollen  + wingmm + (1 |Colony/ Microcol), 
                  family = "nbinom",
                  data=Winged)  
summary (nbinom01) #let's remove wing size
nbinom001<-glmmadmb(totcrith ~ Pollen  +  (1 |Colony/ Microcol), 
                   family = "nbinom",
                   data=Data) #great, no error!
summary(nbinom001) #big effect of pollen

Anova(nbinom001)
summary(nbinom001)
library(bbmle)
AICtab (nbinom01,nbinom001 )
AIC(nbinom01) 
AIC(nbinom001) #they are pretty close
summary(nbinom1)
CountFinal<-nbinom001
#2 methods to test individual terms:
#1. chi-square test
Anova(CountFinal)

#2. likelihood ratio test (alternative to chi-squared)
CountStripped<-update(CountFinal, ~. - Pollen) #makes a reduced model
anova(CountStripped, CountFinal) #compares reduced model to the model with pollen as predictor




#try with a different package for comparison
library(lme4)
?glmer
M0<- glmer.nb(totcrith ~ Pollen + wingmm + dimorph + (1|Colony) + (1|Microcolony), data=Dimorphed)
summary(M0)
Anova(M0)
#drop wing size and dimorphism

M1<-glmer.nb(totcrith ~ Pollen + (1|Colony) + (1|Microcolony) , data=Data)
#compare crossed (above) with explicit nesting structure (below)
M1a<-glmer.nb(totcrith ~ Pollen + (1|Colony/Microcol) , data=Data)

summary(M1)
summary(M1a)
#exactly the same, p-value quite a bit higher(!) than with glmmadmb
summary(M1)
#still get warning
Anova(M1) 

#Compare fits with 2 packages
logLik(M1)
#compare to glmmadmb fit:
logLik(CountFinal) #admb fit much better ("CountFinal")
AICtab (M1, CountFinal)
#glmmadmb better by 80(!) aic units


####PLOTTING
#Use lsmeans to get mean counts
library(lsmeans)
My.Lsmeans <- lsmeans(CountFinal, ~Pollen)
My.Lsmeans
altLsmeans<- lsmeans(M1, ~Pollen)
altLsmeans #somewhat smaller confidence bands
#I am going to go with the glmmadmb model ("CountFinal"), 
#since it did not give any warnings

#BELOW adapted FROM CONROY EXPERIMENT
#######PLOTTING : PARASITE LOAD#############
#plotting marginal means:
library(ggplot2)
#make a data frame from the lsmeans object
Lsmeans.df<-summary(My.Lsmeans)
str(Lsmeans.df)
View(Lsmeans.df)
#label for x axis
pollen.axis<-c( "Buckwheat", "Sunflower")

plot0<- ggplot(Lsmeans.df, aes(x=Pollen, y=lsmean)) + 
  geom_bar(position=position_dodge(), 
           stat="identity", colour="black",
           size=1.5) +
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE),
                size=1.5, width=0.2,  # Width of the error bars
                position=position_dodge(0.9))
plot0 #Solid base plot
plot01<-plot0 + theme_bw() #remove colored background
ylabel<- expression(bold(paste("Parasite load (ln (cells * 0.02 ", mu,L^-1, "))", sep="") )) #check for encoding
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
        panel.background = element_blank()) 
neater
#axis scaling
Zeroed<- neater + coord_cartesian(ylim = c(0, 4.3))
Zeroed
#if we want y axis to start at 0:
#this removes space between x axis and bottom of bar (use for non-log scale)
Justified<- neater +
scale_y_continuous(expand = c(0,0), limits = c(0,4.3))
Justified
#lost axes
J2<-Justified + theme(axis.line.y=element_line(size=2),
                  axis.line.x=element_line(size=2))

#adjust legend
#fig1<-MoveLeg #or "Justified" to start y axis at 0
#fig1<-fig1 +  theme(axis.title.y=element_text(vjust=0.9))
fig1<-J2
#export to PDF
#dev.new() #sometimes this helps with errors
setwd("C:/Users/Evan/Dropbox/Jess Sum15/jess stats/Counts/R_analyses")
pdf("fig1_counts_v2.pdf",height=7, width=5, paper='special') 
fig1
dev.off()

#for eps
setEPS()
postscript("fig1countsv2.eps", horizontal = FALSE, onefile = FALSE, 
           paper = "special", height = 7, width = 5)
fig1
dev.off()



#citation info
setwd("C:/Users/Evan/Dropbox/Jess Sum15/jess stats/Counts/R_analyses")
#install.packages("bibtex")
library(bibtex)
write.bib(entry="ggplot2", file = "ggplot2.bib", append = FALSE, verbose = TRUE)
write.bib(entry="coxme", file = "coxme.bib", append = FALSE, verbose = TRUE)
write.bib(entry="survival", file = "survival.bib", append = FALSE, verbose = TRUE)
write.bib(entry="car", file = "car.bib", append = FALSE, verbose = TRUE)
write.bib(entry="lme4", file = "lme4.bib", append = FALSE, verbose = TRUE)
write.bib(entry="MASS", file = "MASS.bib", append = FALSE, verbose = TRUE)
write.bib(entry="lsmeans", file = "lsmeans.bib", append = FALSE, verbose = TRUE)
write.bib(entry="glmmADMB", file = "glmmADMB.bib", append = FALSE, verbose = TRUE)










