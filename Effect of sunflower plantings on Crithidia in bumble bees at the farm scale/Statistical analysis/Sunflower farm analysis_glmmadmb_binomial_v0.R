#######SCRIPT TO ANALYZE CRITHIDIA FROM WILD BEES########
###Collected from Amherst area farms summer 2015
####PROPOSED EXPLANATORY VARIABLE: sunflower area#

#read data
#Google Drive (epy)
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Analysis.Individual.Bee.and.Farm/Analyses/Farm Sampling")
farmdata<-read.csv("Analysis_1_All.csv")



#JG machine:
farmdata<-read.csv("/Users/rg9403/Google Drive/Project_Sunflower/Analysis.Individual.Bee.and.Farm/DATA/CSV/Farm Sampling/Analysis_1_All.csv",header=TRUE)
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
#replace Count>0 with "1" for binomial analysis
farmdata$Infection<-farmdata$Count
farmdata$Infection[farmdata$Infection>0] <- 1
farmdata$Infection
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

plot(Data$Infection~Data$SunArea)
byarea<-lm(Data$Infection~Data$SunArea)
abline(byarea)
summary(byarea) 
Anova(byarea) #interesting trend there
#but it is not accounting for the non-independence
#of bees from same farm

plot(Data$Infection~Data$Log1P.SunArea.)
bylog<-lm(Data$Infection~Data$Log1P.SunArea.) #steep!
abline(bylog)
summary(bylog) #interesting, quite significant trend
#but does not account for non-independence 
#of bees from same farm

#here let's also scout the effects of collection date
plot(Data$CollectionDate,Data$Infection)
#decreasing crithidia over time
bydate<-lm(Data$Infection~Data$CollectionDate)
abline(bydate)
Anova(bydate) #steep decline over time


#let's run a binomial model
library(lme4)

M0<-glmer(Infection~ Log1P.SunArea. + WingSize + CollectionDate+
        (1|Farm),  family = "binomial",
      nAGQ = 25, data=Data) #convergence failure
summary(M0)

M0Laplace<-glmer(Infection~ Log1P.SunArea. + WingSize + CollectionDate+
            (1|Farm),  family = "binomial",
           data=Data) #convergence failure
noWing<-update(M0, ~. - WingSize)
summary(noWing)
anova(M0, noWing) #first model better
noDate<-update(M0, ~. -CollectionDate)
#ahh, no errors
anova(M0, noDate)
#but first model was better!

summary(M0)

NoLog<-glmer(Infection~SunArea + WingSize + CollectionDate+
               (1|Farm),  family = "binomial",
             data=Data) 
#more warnings!


#see if glmmADMB can handle this
library(glmmADMB)
M0g<-glmmadmb(Infection ~ SunArea + WingSize + CollectionDate+
                (1|Farm),  family = "nbinom",
              zeroInflation = TRUE, data=Data)
#had some trouble, unable to estimate standard errors
summary(M0g)


#somehow Jon got this script to run ... ?without errors??
farmdata$Infection1<-as.factor(farmdata$Infection)
TEST<-glmer(Infection~Log1P.SunArea.  + WingSize + CollectionDate+
              (1|Farm),  family = "binomial",
            data=farmdata) 
summary(TEST)


Testg<-glmmadmb(Infection~Log1P.SunArea.  + WingSize + CollectionDate+
                  (1|Farm),  family = "binomial",
                data=Data) 
summary(Testg)
Anova(Testg)
# Analysis of Deviance Table (Type II tests)
# 
# Response: Infection
# Df  Chisq Pr(>Chisq)   
# Log1P.SunArea.   1 3.3278   0.068119 . 
# WingSize         1 6.2236   0.012606 * 
#   CollectionDate   1 8.2538   0.004067 **
#   Residuals      655    

library(lsmeans)
#make a reference grid
#see worked example at bottom of
#http://www.ats.ucla.edu/stat/r/dae/zinbreg.htm
summary(TEST)
#example:
#newdata1 <- expand.grid(0:3, factor(0:1), 1:4)
#this creates a grid of predictor variables
#means "create grid of all combos with
#first column 0 to 3, second column 0 and 1, third column 1 to 4... 
#colnames(newdata1) <- c("child", "camper", "persons") 
#gives each column a name
#newdata1$phat <- predict(m1, newdata1) #make col. for predicted value
#computes a column for response variable based on predictors
areas<-seq(0,max(Data$Log1P.SunArea.), 0.01)
dates<-mean(Data$CollectionDate)
Wing<-mean(Data$WingSize, na.rm=TRUE)
newdata<-expand.grid(areas, dates, Wing)

colnames(newdata) <- c("Log1P.SunArea.", "CollectionDate", "WingSize") 
newdata

?predict


output<-predict(Testg, newdata, se.fit=TRUE)


#cool, get predictions w/ standard errors!
newdata$phat<-output$fit
newdata$se<-output$se.fit
newdata$plus<-newdata$phat + newdata$se
newdata$minus<-newdata$phat - newdata$se
newdata$CIup<-newdata$phat + 1.96*newdata$se
newdata$CIlo<-newdata$phat - 1.96*newdata$se
newdata$Area=newdata$Log1P.SunArea.

#plot on logit scale (odds ratios)
library(ggplot2)
p<- ggplot(newdata, aes(x=Area, y=phat)) + 
  geom_line(size=2) +
  theme_classic(base_size=24) + 
  theme(axis.line.y=element_line(size=2, color="black"),
        axis.line.x=element_line(size=2, color="black"))
p
onebandSE<-p+geom_ribbon(aes(ymin = minus, ymax = plus), alpha = .1)
onebandSE
#some pretty big bands!



#back transform to scale of response
#install.packages("gtools")
library(gtools)
newdata$phatRAW<-inv.logit(newdata$phat)
newdata$plusSERAW<-inv.logit(newdata$plus)
newdata$minusSERAW<-inv.logit(newdata$minus)
newdata$CIupRAW<-inv.logit(newdata$CIup)
newdata$CIloRAW<-inv.logit(newdata$CIlo)

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
scattered_v0<-p0 + geom_point(data=Data,
                                  aes(x=Log1P.SunArea., y=Infection))
xlabel<-
  expression(bold('Log'[10]* " Sunflower area (acres)"))
labeled<-scattered_v0 + ylab("Infection presence") + xlab(xlabel) + 
  theme(text=element_text(face="bold"))
labeled

#jpeg("farm_binomial.jpg")
pdf("farm_binomial.pdf", height=7, width=7, paper="special")
labeled
dev.off()






#add confidence bands
#this was not all that informative----
oneband<-p0+geom_ribbon(aes(ymin = CIloRAW, ymax = CIupRAW), alpha = .1)
onebandSE<-p0+geom_ribbon(aes(ymin = minusSERAW, ymax = plusSERAW), alpha = .1)
oneband  #considerable variation!
onebandSE #still all over the place
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
