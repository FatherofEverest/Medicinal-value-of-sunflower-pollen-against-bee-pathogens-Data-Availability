####analysis of POLLEN EXPERIMENT_2_JC Roulston Crithidia Strain ####### 
###individual bees.... JONATHAN GIACOMINI.... 3 POLLEN (Sun + Buck + KMix) DIETS EFFECTS ON CRITHIDIA COUNTS#
### All bees were inoculated with Crithidia bombi on day 1 and day 2 (double dose)
### Bees were fed Kmix for first week and either switched to Buck or Sun, or remained on Kmix, for an additional week

data<-read.csv("/Users/rg9403/Dropbox/Sunflower Pollen Experiments/Sunflower Pollen Exp 2_Round2.csv",header=TRUE)
str(data) 

###Need to re-install glmmADMB, Not sure why, maybe there was an update required
install.packages("R2admb")
install.packages("glmmADMB", 
                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
                         getOption("repos")),
                 type="source")
library(glmmADMB)


############## Negative Binomial Fit ##################
Mod1<-glmmadmb(GutCount ~ Treatment  + WingSize +(1|ColonyID), 
                  family = "nbinom",
                  data=data)
summary(Mod1)
library(car)
Anova(Mod1)
# Analysis of Deviance Table (Type II tests)
# 
# Response: GutCount
# Df   Chisq Pr(>Chisq)    
# Treatment  2 23.9061  6.439e-06 ***
#   WingSize   1  1.6613     0.1974    
# Residuals 43   

### WingSize not sig (P = 0.1974), Let's drop it fromt he model
Mod2<-update(Mod1, ~. - WingSize)
summary(Mod2)
Anova(Mod2)
# Analysis of Deviance Table (Type II tests)
# 
# Response: GutCount
#           Df  Chisq   Pr(>Chisq)    
# Treatment  2 22.086    1.6e-05 ***
#   Residuals 44 

AIC(Mod1) #470.988
AIC(Mod2) #470.604 Not much difference 

#post hoc test 
library(multcomp)
summary(glht(Mod2,linfct=mcp(Treatment="Tukey")))
# Linear Hypotheses:
#                     Estimate Std. Error z value Pr(>|z|)    
# KMIX - BUCK == 0  -0.3631     0.5324  -0.682   0.7634    
# SUN - BUCK == 0   -2.4561     0.5291  -4.642   <0.001 ***
#   SUN - KMIX == 0   -2.0930     0.7545  -2.774   0.0142 * 

######## PLOT ########
library(lsmeans)
My.Lsmeans <- lsmeans(Mod2, ~Treatment) 
My.Lsmeans
# Treatment   lsmean  SE       df asymp.LCL asymp.UCL
# BUCK      4.826868 0.3763100 NA  4.089314  5.564422
# KMIX      4.463766 0.6494017 NA  3.190962  5.736570
# SUN       2.370775 0.6538854 NA  1.089183  3.652367

library(ggplot2)
Lsmeans.df<-summary(My.Lsmeans)
str(Lsmeans.df)
View(Lsmeans.df)
#label for x axis
pollen.axis<-c( "Buckwheat","K-Mix", "Sunflower")

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

