#######SCRIPT TO ANALYZE CRITHIDIA FROM WILD BEES########
###Collected from Amherst area farms summer 2015
####PROPOSED EXPLANATORY VARIABLE: sunflower area#
rm(list=ls())
#read data
#Google Drive (epy)
setwd("C:/Users/Evan/Google Drive/ The Sunflower Pollen Project/Analysis.Individual.Bee.and.Farm/Analyses/Farm Sampling")
farmdata<-read.csv("Analysis_1_All.csv")

#JG machine:
farmdata<-read.csv("/Users/rg9403/Desktop/sunflowersamplepatrickdata/SUNFLOWER_FARM_SAMPLING/Sunflower Sampling Data_CSVs/Analysis_1_All.csv",header=TRUE)
str(farmdata)
# 'data.frame':	667 obs. of  13 variables:
#   $ CollectionDate : int  : Julian date of collection
# $ Farm           : Factor w/ 22 levels : denotes farm where bees collected
# $ BeeID  : unique number for each bee      
# $ WingSize       : num  : size of marginal cell of right forewing
# $ FlowerID       : Factor w/ 18 levels "BAS","BRG","BSUN",..: 2 3 3 3 3 3 7 9 9 9 ...
# $ Count          : int  : count of motile Crithidia cells in 0.02uL of gut supernatant
# $ Log1P.Count.   : log10 of "Count"?
# $ FlowerID2 ????means what???     : Factor w/ 2 levels "NONSUN","SUN": 1 2 2 2 2 2 1 1 1 1 ...
# $ FlowerID3 ???means what?     : Factor w/ 3 levels ".","HELI","NONHELI": 3 2 2 2 2 2 3 3 3 3 ...
# $ FarmSize  ???units??     : Factor w/ 22 levels "112,137.41","131,041.50",..: 10 10 10 10 10 10 10 10 10 10 ...
# $ Log1P.FarmSize. ???units???: num  4.52 4.52 4.52 4.52 4.52 4.52 4.52 4.52 4.52 4.52 ...
# $ SunArea        : num  67.9 67.9 67.9 67.9 67.9 ...
# $ Log1P.SunArea. : 
###     EPY COMMENT-- I HAVE A PROBLEM WITH "Log1P.SunArea"
###There are zeroes for "SunArea"
#### Check correlations between explanatory variables, especially collectiondate.####
d<- data.frame(farmdata$CollectionDate,farmdata$Log1P.SunArea.,farmdata$WingSize)
cor(d,use="complete.obs")
#                          farmdata.CollectionDate farmdata.Log1P.SunArea. farmdata.WingSize
# farmdata.CollectionDate               1.0000000              0.02159180        0.18897770
# farmdata.Log1P.SunArea.               0.0215918              1.00000000       -0.08091181
# farmdata.WingSize                     0.1889777             -0.08091181        1.00000000
#none of the variables strongly correlated

#farmdata$CollectionDate<-factor(farmdata$CollectionDate)
#EPY note-- collection date I will keep as numeric covariate....
#to account for increasing/decreasing infection across the sampling period
farmdata$Count<-as.numeric(farmdata$Count)

Data<-na.omit(farmdata) #in order for glmmadmb's to run

Data1<-Data
###Scout some distributions, use Julia Pilowky script
#http://www.juliapilowsky.com/mixedmodels/
#normal
library(MASS)
library(car)
qqp(Data1$Count, "norm") #overdispersed
# lnorm means lognormal
qqp(Data1$Count, "lnorm") 
#for this we would have to transform zeroes to ones

#Poisson:
poisson <- fitdistr(Data1$Count, "Poisson")
qqp(Data1$Count, "pois", poisson$estimate)
#not good, way over-dispersed

#negative binomial
nbinom <- fitdistr(Data1$Count, "Negative Binomial")
qqp(Data1$Count, "nbinom", size = nbinom$estimate[[1]],
    mu = nbinom$estimate[[2]])
#that's pretty terrible, under-dispersed
nbinom <- fitdistr(1+Data1$Count, "Negative Binomial")
qqp(1+Data1$Count, "nbinom", size = nbinom$estimate[[1]],
    mu = nbinom$estimate[[2]]) #worked better

plot(Data$Count~Data$SunArea)
byarea<-lm(Data$Count~Data$SunArea)
abline(byarea)
summary(byarea) 
Anova(byarea) #interesting trend there
#but it is not accounting for the non-independence
#of bees from same farm

plot(Data$Count~Data$Log1P.SunArea.)
bylog<-lm(Data$Count~Data$Log1P.SunArea.) #steep!
abline(bylog)
summary(bylog) #interesting, quite significant trend
#but does not account for non-independence 
#of bees from same farm

#here let's also scout the effects of collection date
plot(Data$CollectionDate,Data$Count)
#decreasing crithidia over time
bydate<-lm(Data$Count~Data$CollectionDate)
abline(bydate)
Anova(bydate)


################################ ZERO-INFLATED NEGATIVE BINOMIAL ####################################
library(glmmADMB)
#epy note 7 march: using raw area, not log scale here
#also include julian date as covariate
M0<-glmmadmb(Count ~ SunArea + WingSize + CollectionDate+
               (1|Farm),  family = "nbinom",
             zeroInflation = TRUE, data=Data)
summary(M0) #highly significant effect of collection date

#note the highly significant zero inflation (I think):
#Zero-inflation: 0.48902  (std. err.:  0.025489 )

#co-authors are not going to go for that...
#beause we are invoking sunflower as creator of the zeroes
M1<-glmmadmb(Count ~ SunArea + WingSize + CollectionDate+
               (1|Farm),  family = "nbinom",
             zeroInflation = FALSE, data=Data)
summary(M1) #highly significant effect of collection date
library(bbmle)
AICtab(M0, M1) #M0 is a better fit, 
#probably because we have too many zeroes for the regular model
#but we want to analyze the zeroes and attribute to sunflower

####so I (EPy) am feeling pulled toward the binomial model now (7 march 2016)
#this was the original model coded by JG
#EPY don't understand why Log(0) area is computed as 0, 
#should be negative infinity?
Log0<- glmmadmb(Count ~ Log1P.SunArea.  + WingSize + 
                   CollectionDate +
                   (1|Farm), 
                 family = "nbinom",zeroInflation = FALSE, data=Data) 
summary(Log0)
logLik(Log0)
logLik(M0)
logLik(M1)
AICtab(M0, M1, Log0) #the zero-inflated model 
#is still working better for this distribution, but 
#perhaps not for our needs

#Let's drop wing size p=0.6
L1<-update(Log0, ~. -WingSize)
anova(Log0, L1)
AICtab(Log0, L1) #keep the simpler model
str(Data$WingSize)

#now we can reinstate the bees
#for whom we did not have wing measures
L1plus<-update(L1, data=farmdata)
#looks good-- all predictors are helping us there
summary(L1plus)
Anova(L1plus)

#so this is the model of choice:
#L1plus<-glmmadmb(Count ~ Log1P.SunArea. + CollectionDate +(1|Farm), 
        #family = "nbinom",zeroInflation = FALSE, data=farmdata) 

library(effects)
allEffects(L1plus) #unhappy

library(lsmeans)
#make a reference grid
#see worked example at bottom of
#http://www.ats.ucla.edu/stat/r/dae/zinbreg.htm
summary(L1plus)
#example:
#newdata1 <- expand.grid(0:3, factor(0:1), 1:4)
#this creates a grid of predictor variables
#means "create grid of all combos with
#first column 0 to 3, second column 0 and 1, third column 1 to 4... 
#colnames(newdata1) <- c("child", "camper", "persons") 
#gives each column a name
#newdata1$phat <- predict(m1, newdata1) #make col. for predicted value
#computes a column for response variable based on predictors
areas<-seq(0,max(farmdata$Log1P.SunArea.), 0.1)
dates<-mean(farmdata$CollectionDate)
newdata<-expand.grid(areas, dates)
colnames(newdata) <- c("Log1P.SunArea.", "CollectionDate") 
newdata

library(lsmeans)
?predict
lsmeans(L1plus, ~Log1P.SunArea., 
        at=areas)

output<-predict(L1plus, newdata, se.fit=TRUE)
#cool, get predictions w/ standard errors!
newdata$phat<-output$fit
newdata$se<-output$se.fit
newdata$plus<-newdata$phat + newdata$se
newdata$minus<-newdata$phat - newdata$se
newdata$CIup<-newdata$phat + 1.96*newdata$se
newdata$CIlo<-newdata$phat - 1.96*newdata$se

#back transform to scale of response
newdata$phatRAW<-exp(newdata$phat)
newdata$plusSERAW<-exp(newdata$plus)
newdata$minusSERAW<-exp(newdata$minus)
newdata$CIupRAW<-exp(newdata$CIup)
newdata$CIloRAW<-exp(newdata$CIlo)

View(newdata)

#rename column for convenience
newdata$Area=newdata$Log1P.SunArea.

library(ggplot2)
p0<- ggplot(newdata, aes(x=Area, y=phatRAW)) + 
  geom_line(size=2) +
 theme_classic(base_size=24) + 
  theme(axis.line.y=element_line(size=2, color="black"),
        axis.line.x=element_line(size=2, color="black"))
p0
#add confidence bands
oneband<-p0+geom_ribbon(aes(ymin = CIloRAW, ymax = CIupRAW), alpha = .1)
onebandSE<-p0+geom_ribbon(aes(ymin = minusSERAW, ymax = plusSERAW), alpha = .1)
oneband  #considerable variation!
onebandSE #a little more manageable
#now add the raw data
scattered<-onebandSE + geom_point(data=farmdata,
                aes(x=Log1P.SunArea., y=Count))
scattered
ylabel<-expression(bold(paste("Parasite load (cells * 0.02 ",
                        mu, L^-1, "))", sep="") )) #check for encoding
xlabel<-
  expression(bold('Log'[10]* " Sunflower area (acres)"))
scattered + xlab(xlabel)
 # xlab=expression('hi'[5]*'there'[6]^8*'you'['down here']*'and'^'up'*'there'))
Fig1<-scattered + xlab(xlabel) +
  #ylab("Parasite load")
  ylab(ylabel) + theme(text=element_text(face="bold"))
Fig1

jpeg("Farmfig.jpg")
#pdf("Farmfig.pdf", height=7, width=7, paper="special")
Fig1
dev.off()
