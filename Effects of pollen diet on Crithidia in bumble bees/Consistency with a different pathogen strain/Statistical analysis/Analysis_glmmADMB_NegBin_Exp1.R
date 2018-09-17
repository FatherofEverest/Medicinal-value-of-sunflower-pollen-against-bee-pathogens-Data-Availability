####analysis of POLLEN EXPERIMENT_1_JC Roulston Crithidia Strain ####### 
###individual bees.... JONATHAN GIACOMINI.... 2 POLLEN (Sun + Buck ) DIETS EFFECTS ON CRITHIDIA COUNTS#

data1<-read.csv("/Users/rg9403/Dropbox/Sunflower Pollen Experiments/Exp1/Sunflower Pollen Exp 1_Analysis 1 data.csv",header=TRUE)
str(data1) 

library(glmmADMB)
############## Negative Binomial Fit ##################
data1$Date<-factor(data1$Date) # Change Date to a factor
Data1<-na.omit(data1) #The following models don't work unless I omit the NAs

Fit1<-glmmadmb(Count ~ Treatment  + WingSize +(1|ColonyID)+(1|Date), 
               family = "nbinom",
               data=Data1)
summary(Fit1)
library(car)
Anova(Fit1)
# Analysis of Deviance Table (Type II tests)
# 
# Response: Count
# Df   Chisq Pr(>Chisq)    
# Treatment   1 30.1764  3.945e-08 ***
#   WingSize    1  9.1121   0.002539 ** 
#   Residuals 141

# We can knockout Date and colony ID and then compare models
Fit2<-update(Mod1, ~. - Date)
summary(Fit2)
Anova(Fit2)

Fit3<-update(Mod1, ~. - ColonyID)
summary(Fit3)
Anova(Fit3)


AIC(Fit1) #459.914 #### Best fit...Let's keep Both Colony ID and Date as random effects
AIC(Fit2) #470.988
AIC(Fit3) #470.988

######## PLOT ######## Not so good - the Lsmean for sun is negative 
library(lsmeans)
My.Lsmeans2 <- lsmeans(Fit1, ~Treatment) 
My.Lsmeans2
library(ggplot2)
Lsmeans.df<-summary(My.Lsmeans2)
str(Lsmeans.df)
View(Lsmeans.df)
#label for x axis
pollen.axis<-c( "Buckwheat", "Sunflower")

plot0<- ggplot(Lsmeans.df, aes(x=Treatment, y=lsmean)) + 
  geom_bar(position=position_dodge(), stat="identity", colour="black" ) +
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE),
                size=1, width=0.2,  # Width of the error bars
                position=position_dodge(0.9))
plot0 # Looks ok, let's refine it now
plot01<-plot0 + theme_bw() #remove colored background
ylabel<- expression(bold(paste("Parasite load (ln (cells * 0.02 ", L^-1, "))", sep="") )) #check for encoding
plot02<-plot01 + ylab(ylabel) #y label
plot02
#relabel x axis:
plot05<-plot02 + xlab("Pollen type") + # x axis label
  scale_x_discrete(breaks=levels(Lsmeans.df$Treatment) , labels=pollen.axis)
plot05
#Larger font:
plot06<-plot05 + theme(text=element_text (size=25, face="bold") )
plot06
#Remove gridlines
neater<- plot06 +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.line = element_line(size = 1.5, colour = "black"))
neater
#axis scaling
Zeroed<- neater + coord_cartesian(ylim = c(0, 6.5))
Zeroed
#if we want y axis to start at 0:
#this removes space between x axis and bottom of bar (use for non-log scale)
Justified<- neater +
  scale_y_continuous(expand = c(0,0), limits = c(0,6.5))
Justified
fig1<-Justified
fig1



