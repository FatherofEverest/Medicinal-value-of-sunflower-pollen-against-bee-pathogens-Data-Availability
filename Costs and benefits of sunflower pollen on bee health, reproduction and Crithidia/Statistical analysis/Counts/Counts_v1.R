
#Analysis of Crithidia counts in microcolonies fed different pollens
#Jessica Leslie honors thesis 2015

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

#load packages
#run "install.packages("packagename") if you need to install-- remember the quotation marks!
library("glmmADMB") #for mixed model (glmmadmb)
library(lme4) #alternative mixed model, function "glmer"
library("ggplot2") #for plotting
library("car") #chi-square tests (Anova)
library("lsmeans")


#set working directory
setwd("C:/Users/Evan/Dropbox/Jess Sum15/jess stats/Counts/R_analyses")
Data<-read.csv("countsforJMP.csv")
View(Data)
#install.packages("R2admb") #run if need to install, takes some time esp. next line
#install.packages("glmmADMB", 
repos=c("http://glmmadmb.r-forge.r-project.org/repos",
        getOption("repos")),
type="source")

#check dataset
str(Data)
Data$Microcol<-factor(Data$Microcol)
#eliminate bees without counts
Counted<-subset(Data, !totcrith=="NA")
#strip microcolonies without dimorphism
Counted1<-subset(Counted, !dimorph=="NA")
Counted2<-subset(Counted1, !wingmm=="NA")

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
                    data=Counted2)  #throwing errors with poisson fit, no good
#now negative binomial
#Start with full nesting
nbinom0<-glmmadmb(totcrith ~ Pollen  +  wingmm + dimorph + (1 |colony/ Microcol/BeeID), 
                  family = "nbinom",
                  data=Counted2)  
Anova(nbinom0)
#eliminate "dimorph" 
nbinom1<-update(nbinom0, ~. - dimorph) #drop dimorphism
anova(nbinom1, nbinom0) #no difference
#now drop wingmm
nbinom2<-update(nbinom1, ~. - wingmm)
anova(nbinom2, nbinom1) #p=0.1456
Anova(nbinom2) #p=0.1458
library(  "bbmle")
AICtab (nbinom2, nbinom1) #minimal difference, so elminate wing size
summary(nbinom2)
#now we can add back bees missing wing measurements by using data=Counted:
CountFinal<-glmmadmb(formula = totcrith ~ Pollen + (1 | colony/Microcol/BeeID), 
         data = Counted, family = "nbinom")
#2 methods to test individual terms:
#1. chi-square test
Anova(CountFinal)
#Analysis of Deviance Table (Type II tests)
#Response: totcrith
#Df Chisq Pr(>Chisq)    
#Pollen  1 41.76  1.032e-10 ***
#2. likelihood ratio test
CountStripped<-update(CountFinal, ~. - Pollen) #makes a reduced model
anova(CountStripped, CountFinal) #compares reduced model to the model with pollen as predictor
#Analysis of Deviance Table
#Model 1: totcrith ~ 1
#Model 2: totcrith ~ Pollen
#NoPar LogLik Df Deviance  Pr(>Chi)    
#1     5 -570.2                          
#2     6 -554.4  1     31.6 1.894e-08 ***

nbinom1<-update(nbinom, ~. - dimorph) #drop dimorphism
Anova(nbinom1)
nbinom2<-update(nbinom1, ~. - wingmm)
Anova(nbinom2)
nbinom2
#try zero inflation:
zinfl<- glmmadmb(totcrith ~ Pollen + (1 | Microcol), data=Counted2,
                 family = "nbinom",zeroInflation = TRUE )
#oops-- got error. Stick with regular negative binomial for now

#try nesting of random effect
summary(nbinom2)
#I used data=Counted to add back the bees missing wing size or dimorphism measures
nested<- glmmadmb(formula = totcrith ~ Pollen + (1 |colony/ Microcol), data = Counted, 
         family = "nbinom")
summary(nested)


#try to add one more level of nesting for BeeID
nested2<- glmmadmb(formula = totcrith ~ Pollen + (1 |colony/ Microcol/BeeID), data = Counted, 
                   family = "nbinom")
summary(nested2)



#try with a different package

library(lme4)
?glmer
M1<-glmer.nb(totcrith ~ Pollen + (1|Microcol), data=Counted)
#still get warning
Anova(M1)


#use "nbinom2" as final model-- no warnings from that model
Anova(nbinom2)


glmer.nb(totcrith ~ Pollen  + #fixed
           daysinfect +dimorph + 1|colony/Microcol/BeeID, 
         zeroInflation = TRUE, 
         data=Counted)



fitPoiss<- glmmadmb(totcrith ~ pollen  + #fixed
                      daysinfect + dimorph +
                      1( colony|Microcol), 
                    zeroInflation = TRUE, family = "poisson",
                    data-Data) 
ezfit<-glmmadmb(totcrith~pollen + )
?glmmadmb
