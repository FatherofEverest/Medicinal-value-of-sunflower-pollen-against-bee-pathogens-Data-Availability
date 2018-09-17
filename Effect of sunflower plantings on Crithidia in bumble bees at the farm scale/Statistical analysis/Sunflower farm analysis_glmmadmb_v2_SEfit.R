#######SCRIPT TO ANALYZE CRITHIDIA FROM WILD BEES########
###Collected from Amherst area farms summer 2015
####PROPOSED EXPLANATORY VARIABLE: sunflower area#
rm(list=ls())
#read data
#Google Drive (epy)
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Analysis.Individual.Bee.and.Farm/Analyses/Farm Sampling")
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

#simple boxplot
library(car)
Boxplot(Count~Log1P.SunArea., data = Data1)


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
plot(M0$residuals)
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
Log0<- glmmadmb(Count ~ Log1P.SunArea.  + WingSize + 
                   CollectionDate +
                   (1|Farm), 
                 family = "nbinom",zeroInflation = FALSE, data=Data) 
summary(Log0)
#Compare to zero-inflated:
Log.zi<-glmmadmb(Count ~ Log1P.SunArea.  + WingSize + 
                   CollectionDate +
                   (1|Farm), 
                 family = "nbinom",zeroInflation = TRUE, data=Data) 
library(pscl)
vuong(Log0, Log.zi) #error
library(bbmle)
AICtab(Log0, Log.zi) #this is quite a big improvement with zero inflation
#Tough call whether there are 2 processes at work or not
#In some sense, every bee ought to have been exposed to some Crithidia 
  #by the time it starts to forage
#However, one could envision that some bees have not had sufficient time or exposure
   #to develop detectable infection
plot(Log0$residuals)
plot(Log0$fitted, Log0$residuals)
qqp(Log0$residuals) #oof

plot(Log.zi$fitted, Log.zi$residuals)
plot(Log.zi$residuals)
qqp(Log.zi$residuals) #this model also veers off the line at high values

summary(Log0)
summary(Log.zi) #Zero-inflation: 0.4862  (std. err.:  0.026077 )


logLik(Log0)
logLik(M0)
logLik(M1)
AICtab(M0, M1, Log0) #the zero-inflated model 
#is still working better for this distribution, but 
#perhaps not for our needs. Can't decide. For now use the non-zi model for consistency with lab expts

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

#############PLOTTING##############
library(lsmeans)
MyModel<-L1plus
areas<-seq(-0.2,max(0.1+farmdata$Log1P.SunArea.), 0.1)
####scale of response############
ls2plot<- summary(lsmeans(MyModel, ~Log1P.SunArea., 
                          at=list(Log1P.SunArea.=areas),
                          type="response"))
#View(ls2plot)
########margina means on scale of linear predictor##########
margmeans<-summary(lsmeans(MyModel, ~Log1P.SunArea., 
                           at=list(Log1P.SunArea.=areas)))
#Manually exponentiate the (mean+-SE) from scale of linear predictor
margmeans$btmean<-exp(margmeans$lsmean)
margmeans$rawhi<-margmeans$lsmean + margmeans$SE
margmeans$rawlo<-margmeans$lsmean - margmeans$SE
margmeans$bthi<-exp(margmeans$rawhi)
margmeans$btlo<-exp(margmeans$rawlo)                        
margmeans$Area=margmeans$Log1P.SunArea.


#rename column for convenience
ls2plot$Area=ls2plot$Log1P.SunArea.

library(ggplot2)
########OPTION 1: SYMMETRIC ERROR BARS WITH "SE" from lsmeans type="response"####
##note this will give part of confidence bands below 0 ??!!######
p0<- ggplot(ls2plot, aes(x=Area, y=response)) + 
  geom_line(size=2) + theme_classic(base_size=24) + 
  theme(axis.line.y=element_line(size=2, color="black"),
        axis.line.x=element_line(size=2, color="black"))
#add confidence bands
p01<-p0+geom_ribbon(aes(ymin=response - SE, ymax=response + SE), alpha=0.5) #pretty!
p01



########OPTION 2 (preferred):#################
########Asymmetric error bands from manually exponentiated mean +-SE
p0<- ggplot(margmeans, aes(x=Area, y=btmean)) + 
  geom_line(size=2) + theme_classic(base_size=24) + 
  theme(axis.line.y=element_line(size=2, color="black"),
        axis.line.x=element_line(size=2, color="black"))
  #add confidence bands
  p01<-p0+geom_ribbon(aes(ymin=btlo, ymax=bthi), alpha=0.5) #pretty!
p01
ylab<-expression(bold(Parasite~load~("cells*"~0.02~mu*L^-1)))
xlabel<-  expression(bold(Log[10]~(Sunflower~area~(m^2))))
p1<-p01 +  theme_classic(base_size = 25) +
  xlab(xlabel)+ylab(ylab)
p2<-p1 +  theme(text=element_text(face="bold"), 
        line=element_line(size=2)) +
  theme(axis.line.x=element_line(size=2),
        axis.line.y=element_line(size=2))
p3<- p2 + scale_x_continuous(expand = c(0,0))#removes gap
p3

#now add the raw data
scattered<-p2 + geom_point(data=farmdata,
                aes(x=Log1P.SunArea., y=Count))
scattered

#jpeg("Farmfig_v3_manualSE.jpg")
#pdf("Farmfig_v3_manualSE.pdf", height=7, width=7, paper="special")
scattered
#dev.off()

#ggsave("Farmfig_v3_manualSE.eps",height=7, width=7, device=cairo_ps)

################################################################
####update 2017.01
##jg requested plot with tighter y-axis scale
#a log-scale would make the model appear less diffuse
#####################################################################

#take your pick of scaling:

(lnplot.autoscale<-scattered + scale_y_continuous(trans = "log1p", breaks = waiver()))


log(exp(1))

scattered + scale_y_continuous(trans = "log1p", breaks = c(exp(1), exp(2), exp(3), exp(4), exp(5), exp(6)))
b1<- expression(e^1)
b2<- expression (e^2)
b3<- expression (e^3)
b4<- expression (e^4)
b5<- expression (e^5)
b6<- expression (e^6)
(lnplot.ESCALE<-scattered + scale_y_continuous(trans = "log1p", breaks = c(exp(1), exp(2), exp(3), exp(4), exp(5), exp(6)),
                               labels = c(b1, b2, b3, b4, b5, b6)))
#pdf("farmplot.logscale.ELABELS.pdf")
lnplot.ESCALE
 #    dev.off()                          



######
#make plot from scratch using ln-scale:
#####                                   ##############
View(ls2plot)
lnp1<- 
  ggplot(margmeans, aes(x=Area, y=lsmean)) + 
  geom_line(size=2) + theme_classic(base_size=24) + 
  theme(axis.line.y=element_line(size=2, color="black"),
        axis.line.x=element_line(size=2, color="black"))
#add confidence bands
lnp2<-lnp1+geom_ribbon(aes(ymin=lsmean - SE, ymax=lsmean + SE), alpha=0.5) #pretty!
lnp2
ylab<-expression(bold(Parasite~load~(ln("cells*"~0.02~mu*L^-1))))
xlabel<-  expression(bold(Log[10]~(Sunflower~area~(m^2))))

lnp3<- lnp2 + theme_classic(base_size = 25) +
  xlab(xlabel)+ylab(ylab)+
  theme(text=element_text(face="bold", color = "black"), 
                line=element_line(size=2, color = "black"),
              axis.line.x=element_line(size=2, color = "black"),
             axis.line.y=element_line(size=2, color = "black"),
              axis.text = element_text(face="bold", color = "black"))

lnp3

###rawdata###
farmdata$lncount<- log(farmdata$Count)
farmdata$lncount #negative infinity for zeroes
farmdata$lncount.tf<- log1p ( farmdata$Count) ####how much to add before the transformation is totally arbitrary
(lnscat<-lnp3 + geom_point(data=farmdata,
                  aes(x=Log1P.SunArea., y=lncount.tf))+
                  theme(legend.position = "none"))

# pdf("farmplot.logscale.pdf", height = 7, width = 7)
# lnscat
# dev.off()

