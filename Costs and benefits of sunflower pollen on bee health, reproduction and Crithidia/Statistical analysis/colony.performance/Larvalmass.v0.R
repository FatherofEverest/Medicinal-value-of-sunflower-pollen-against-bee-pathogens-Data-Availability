####2016.03#########
#analysis of sunflower microcolony performance

#####Larval production################
####Source data: use summary created by LS Adler
#EPY added back colonies from which new adults emerged -- 
#colonies 3, 11, 23, 25, 33
#coded as "1" for pupal production

rm(list=ls()) #clear memory

library(dplyr)
library(car)
library(nlme)
library(lme4)
library(glmmADMB) #takes some work to install this
library(bbmle)
library(ggplot2)

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

#scout the variable
Nothing<-subset(Perf03, larvwt=NA)
NoL <- Perf03[is.na(Perf03$larvwt),]
NoL$larvwt
length(NoL$larvwt)
xtabs(~Pollen, Perf03)
xtabs(~Pollen, NoL) #21 buckwheat, 14 sunflower
xtabs(~Infection, NoL) #17 uninfected, 18 infected
xtabs(~Colony, NoL) #JML-3 and JML-4 fewer larvae, >50% no larvae


Makers<-subset(Perf03, !larvwt=="NA")
Makers$larvwt
length(Makers$larvwt) #41 colonies w/ larvae

#LSA recommends analyzing larval count and larval mass separately
hist(Perf03$larvnum)
LCOUNTS<-ddply(Perf03,~Infection*Pollen, 
          summarise,numlarvae=mean(larvnum),
          meanmass=mean(larvwt, na.rm=TRUE))
LCOUNTS

###HERE I WILL DEAL WITH THE MICROCOLONIES THAT MADE LARVAE
Makers
length(Makers$Microcolony) #41 microcolonies
Boxplot(meanlarvwt~Infection*Pollen, Makers)
#no way that this meets assumptions of equal variance
TABLE<-ddply(Makers,~Infection*Pollen, 
               summarise,numlarvae=mean(larvnum),
               meanmass=mean(larvwt, na.rm=TRUE),
             samplesize=count(larvnum))
xtabs( ~Pollen, Makers)
xtabs( ~Infection, Makers)


#plenty of zeroes
#check distribution
qqp(Perf03$meanlarvwt, "norm") #OK but does not account for different variance!
#anyway the treatment effect is dumb clear in the boxplot,
#so I will just run the model anyway because we have to analyze somehow


#try lognormal
# lnorm means lognormal
qqp(Perf03$meanlarvwt, "lnorm") #we don't need such broad bands
#log-normal would be a biologically appropriate distribution

Buck<-subset(Perf03, Pollen=="B")
Sun<-subset(Perf03, Pollen=="S")

#try gamma
gamma <- fitdistr(Makers$meanlarvwt, "gamma")
gamma <- fitdistr(Makers$meanlarvwt, "gamma")
qqp(Makers$meanlarvwt, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])
#like this b/c allows variance to increase w/ mean


#sample gamma fit
#(m <- glmer(y ~ x + (1 | id), family = Gamma))

library(glmmADMB)
#first try poisson model:
M1<- glmmadmb(meanlarvwt ~ Infection*Pollen  + (1|Colony), 
                    family = "gamma",
                    data=Makers)

summary(M1)
Anova(M1)
M1.ID<-glmmadmb(meanlarvwt ~ Infection*Pollen  + (1|Colony), 
                family = "gamma",link="identity",
                data=Makers) #error
norfit<-glmmadmb(meanlarvwt ~ Infection*Pollen  + (1|Colony), 
                 family = "gaussian",
                 data=Makers) #error
library(bbmle)
AICtab(M1, norfit) #gamma is better
#compare with glmer
M1a<- glmer(meanlarvwt ~ Infection*Pollen  + (1|Colony), 
              family = "Gamma",
              data=Makers) 
summary(M1a)
Anova(M1a) #same idea as with prior model

#drop interaction term
M2<-update(M1, ~. - Infection:Pollen)
summary(M2)
Anova(M2)
M2a<-update(M1a, ~. - Infection:Pollen)
summary(M2a) 
#good-- big effects of pollen
#need to plot
(library(lsmeans))
margmeans<-summary(lsmeans(M2, ~Pollen))
lsmeans(M2, ~Pollen, type="response")

plot(M2$residuals) #looks good, random scatter, no fanning
#awesome! let's plot this-- this will be fun

#get mean, se on scale of response variable
lsmeans(M2, ~Pollen, type="response")
MYMODEL<-M2
#make a data frame from the lsmeans object
Lsmeans.df<-summary(lsmeans(MYMODEL, ~Pollen, type="response"))
#MUST request the "summary" to get a data frame
Lsmeans.df$response


#computing the means +- SE's by manual exponentiation
margmeans$btmean<-exp(margmeans$lsmean)
margmeans$rawhi<-margmeans$lsmean + margmeans$SE
margmeans$rawlo<-margmeans$lsmean - margmeans$SE
margmeans$bthi<-exp(margmeans$rawhi)
margmeans$btlo<-exp(margmeans$rawlo)



##now the plot
library(ggplot2)
#View(Lsmeans.df)
#label for x axis
pollen.axis<-c( "Buckwheat", "Sunflower")
#############option 1###############################
#this gives symmetric error bars using the "SE" column on response scale
plot0<- ggplot(Lsmeans.df, aes(x=Pollen, y=response)) + 
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", size=2 ) +
  geom_errorbar(aes(ymin=response - SE, ymax=response + SE),
                size=2, width=0.2,  # Width of the error bars
                position=position_dodge(0.9))+
  theme_classic(base_size=24)

#############Option 2 (preferred)###################

plot0<- ggplot(margmeans, aes(x=Pollen, y=btmean)) + 
  geom_bar(position=position_dodge(), stat="identity", colour="black", size=2 ) + #color is for OUTLINES
  geom_errorbar(aes(ymin=btlo, ymax=bthi), #adds error bars at 2x thickness
                size=2, width=0.2,  # Width of the error bars
                position=position_dodge(0.9))+ #useful command for line charts to offset error bars, vestigial here?
  theme_classic(base_size=25) #this theme has no background, default base_size is 12 point font

#where are my axes?
p1<-plot0 + ylab("Mean larval mass (g)") + xlab("Pollen")+
   theme (axis.line.y=element_line(color="black", size=2),
           axis.line.x=element_line(color="black", size=2)) +
  theme(text=element_text(face="bold"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.08)) 
p1
p2<-p1+ scale_x_discrete(labels=c("Buckwheat", "Sunflower"))
p2
fig1<-p2

#export to PDF
#dev.new() #sometimes this helps with errors
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Microcolony/jess stats/colony.performance")
pdf("LARVALMASS.v0.manual.SE.pdf",height=6, width=5, paper='special') 
fig1
dev.off()

#for eps
setEPS()
postscript("LARVALMASS.v0.manual.SE.eps", horizontal = FALSE, onefile = FALSE, 
           paper = "special", height = 6, width = 5)
fig1
dev.off()

#just for fun, try another method to reverse transform
#install.packages('RVAideMemoire') #enormous
library(RVAideMemoire)
LSM<-lsmeans(MYMODEL, ~Pollen)
LSM
back.lsmeans(lsm=LSM, transform = c("log") )
exp(-4.47) #not matching the manual calculation...
#seems to be 1 unit lower than the correct estimate
LS1<-summary(LSM)
exp(LS1$lsmean) #differs quite a bit from "back.lsmeans"
#the "lsmeans" package estimates make more sense
#install.packages("DoByMeans") #oops, not available
