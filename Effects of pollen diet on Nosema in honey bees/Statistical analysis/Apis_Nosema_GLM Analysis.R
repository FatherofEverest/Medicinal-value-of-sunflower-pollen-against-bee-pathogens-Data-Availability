############################################################################################################################################################
# Honey Bee - Nosema Exp Analysis
# Performed by JJG - May 25, 2017

#Groups (BeeCupID) of 50 Honey bees (Apis mellifera) were inoculated with Nosema spores and provided either
#sunflower pollen, buckwheat pollen, or no pollen (as a control diet).  Groups of 5 bees and 10 bees were sacraficed on days
#10 and 15, repectively, to measure Nosema infection intensity. 

#A simple GLM should do the trick.  There are no random effects to include in this analysis.  
#Treatment = predictor
#Average.Day15.Spores.ml = Response
#Bonus responses include Day.10.Spores.ml ### See bottom of file 
############################################################################################################################################################
#  Final Model(Day15) = FitNB 
#  Final Model(Day10) = FitNB2  
############################################################################################################################################################
data1<-read.csv("/Users/rg9403/Google Drive/Sunflower Pollen Experiments/Irwin Lab_Sunflower Pollen Experiments/Exp 6 _ Apis mellifera_Nosema_Sun_Buck_NoPollen_June 2016/Honey Bee Nosema Exp_June 2016 (2).csv", header=TRUE)
str(data1)
head(data1)
#   BeeCupID  Treatment  NosemaDay10.Count  Day.10.Spores.mL   Average.Day.15.Count    Average.Day15.Spores.mL
#         1         S               809         40450000                 1491                74550000
#         2         B              1578         78900000                 2035               101750000
#         3        NP               647         32350000                  619                30950000
#         4         S              1018         50900000                 1863                93150000
#         5         B              2048        102400000                 1945                97250000
#         6        NP               774         38700000                  639                31950000

####Renaming to make things no easier
data1$Day10Count<- data1$NosemaDay10.Count
data1$Day15Count<- data1$Average.Day.15.Count
data1$Day10Spores<- data1$Day.10.Spores.mL
data1$Day15Spores<- data1$Average.Day15.Spores.mL
View(data1)

##############################################################################
# Distribution - Response = Day15Spores
##############################################################################
library(car)
library(MASS)

hist(data1$Day15Spores)

qqp(data1$Day15Spores, "norm")
### Not too shabby
qqp(data1$Day15Spores, "lnorm")
### normal fit looks better 

poisson <- fitdistr(data1$Day15Spores, "Poisson")
qqp(data1$Day15Spores, "pois", poisson$estimate)
# Looks nearly the same as the Normal fit


nbinomDist <- fitdistr(data1$Day15Spores, "Negative Binomial")
qqp(data1$Day15Spores, "nbinom", size = nbinom$estimate[[1]],mu = nbinom$estimate[[2]])
### ERRORS: Error in solve.default(res$hessian) : 
###Lapack routine dgesv: system is exactly singular: U[2,2] = 0

##############################################################################
# Nosema Infection intensity models - GLMs
##############################################################################



Fit<- glm(Day15Spores~Treatment, family = gaussian, data = data1)
summary(Fit) #overdispersion because residual deviance much larger than df

FitPoisson<- glm(Day15Spores~Treatment, family = poisson(), data = data1)
summary(FitPoisson) #overdispersion because residual deviance much larger than df

FitQPoisson<- glm(Day15Spores~Treatment, family = quasipoisson(), data = data1)
summary(FitQPoisson) # Again, overdispersion seems to be an issue

# Let's see if any of the models better than the others
anova(Fit,FitPoisson,FitQPoisson, test="Chisq") # Nope



#Ok let's try to wrangle in that overdispersion with a negative binomial model
FitNB<-glm.nb(Day15Spores~Treatment, data = data1)
summary(FitNB) ### Much better. Dispersion seems to be wrangled in now!!
# Call:
#   glm.nb(formula = Day15Spores ~ Treatment, data = data1, init.theta = 32.75755653, 
#          link = log)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.0330  -0.8160  -0.1340   0.7444   2.0697  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) 18.43977    0.05268 350.032  < 2e-16 ***
#   TreatmentNP -1.50116    0.07293 -20.583  < 2e-16 ***
#   TreatmentS  -0.26052    0.07450  -3.497 0.000471 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(32.7576) family taken to be 1)
# 
# Null deviance: 447.726  on 33  degrees of freedom
# Residual deviance:  34.173  on 31  degrees of freedom
# AIC: 1197.3
# 
# Number of Fisher Scoring iterations: 1
# Theta:  32.76 
# Std. Err.:  7.90 
# 2 x log-likelihood:  -1189.305 

Anova(FitNB)
# Analysis of Deviance Table (Type II tests)
# Response: Day15Spores
# LR           Chisq Df Pr(>Chisq)    
# Treatment   413.55  2  < 2.2e-16 ***

summary(glht(FitNB,linfct=mcp(Treatment="Tukey")))

# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glm.nb(formula = Day15Spores ~ Treatment, data = data1, init.theta = 32.75755653, link = log)
# 
# Linear Hypotheses:
#             Estimate  Std. Error   z value     Pr(>|z|)    
# NP - B ==  -1.50116    0.07293     -20.583     < 1e-04 ***
# S -  B ==  -0.26052    0.07450     -3.497      0.00138 ** 
# S - NP ==   1.24064    0.07293      17.011     < 1e-04 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)

Lsmeans.df<-as.data.frame(summary(lsmeans(FitNB, ~Treatment)))
Lsmeans.df$btmean<-exp(Lsmeans.df$lsmean)
Lsmeans.df$btlo<-exp(Lsmeans.df$lsmean - Lsmeans.df$SE) #lower bound
Lsmeans.df$bthi<-exp(Lsmeans.df$lsmean + Lsmeans.df$SE) #upper bound
View(Lsmeans.df)
Lsmeans.df
# Treatment   lsmean         SE df asymp.LCL asymp.UCL    btmean     btlo      bthi
# 1         B 18.43977 0.05268027 NA  18.33652  18.54302 101927273 96696700 107440780
# 2        NP 16.93861 0.05043755 NA  16.83975  17.03747  22716667 21599309  23891827
# 3         S 18.17925 0.05268027 NA  18.07599  18.28250  78550000 74519072  82798971

pollen.axis<-c( "Buckwheat","No-Pollen","Sunflower")
plot0<- ggplot(Lsmeans.df, aes(x=Treatment, y=btmean)) + 
  geom_bar(position=position_dodge(), stat="identity", colour="black" ) +
  geom_errorbar(aes(ymin=btlo, ymax=bthi),
                size=1, width=0.2,  # Width of the error bars
                position=position_dodge(0.9))
plot0
plot01<-plot0 + theme_bw() #remove colored background
ylabel<- expression(bold(paste("Nosema spores/mL", sep="") )) 
plot02<-plot01 + ylab(ylabel) #y label
plot02
#relabel x axis:
plot05<-plot02 + xlab("Treatment") + # x axis label
  scale_x_discrete(breaks=levels(Lsmeans.df$Treatment) , labels=pollen.axis)
plot05

plot06<-plot05 + theme(text=element_text (size=25, face="bold") )
plot06
##remove gridlines (option), remove border
neater<- plot06 + theme(axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), 
                        panel.border = element_blank(),
                        panel.background = element_blank()) + theme(axis.line = element_line(size = 3, colour = "black"))
neater


##############################################################################
# Distribution - Response = Day10Spores
##############################################################################
library(car)
library(MASS)

hist(data1$Day10Spores)

qqp(data1$Day10Spores, "norm")
### Not too shabby # looks a little better than Day15Spores
qqp(data1$Day10Spores, "lnorm")
### not really neccesary

poisson <- fitdistr(data1$Day10Spores, "Poisson")
qqp(data1$Day10Spores, "pois", poisson$estimate)
# Looks nearly the same as the Normal fit

nbinomDist <- fitdistr(data1$Day10Spores, "Negative Binomial")
qqp(data1$Day10Spores, "nbinom", size = nbinom$estimate[[1]],mu = nbinom$estimate[[2]])
### UGGGGHHHHH  ERRORS: Error in solve.default(res$hessian) : 
###Lapack routine dgesv: system is exactly singular: U[2,2] = 0

Fit2<- glm(Day10Spores~Treatment, family = gaussian, data = data1)
summary(Fit2) #overdispersion because residual deviance much larger than df

FitPoisson2<- glm(Day10Spores~Treatment, family = poisson(), data = data1)
summary(FitPoisson2) #overdispersion because residual deviance much larger than df

FitQPoisson2<- glm(Day10Spores~Treatment, family = quasipoisson(), data = data1)
summary(FitQPoisson2) # Again, overdispersion seems to be an issue

# Let's see if any of the models better than the others
anova(Fit2,FitPoisson2,FitQPoisson2, test="Chisq") # Nope

#Ok let's try to wrangle in that overdispersion with a negative binomial model
FitNB2<-glm.nb(Day10Spores~Treatment, data = data1)
summary(FitNB2)
# Call:
#   glm.nb(formula = Day10Spores ~ Treatment, data = data1, init.theta = 18.476068,  link = log)
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.4425  -0.8282  -0.1547   0.5434   3.5404  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) 18.24925    0.07015 260.163  < 2e-16 ***
#   TreatmentNP -0.72006    0.09711  -7.415 1.22e-13 ***
#   TreatmentS  -0.39692    0.09920  -4.001 6.30e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Dispersion parameter for Negative Binomial(18.4761) family taken to be 1)
# Null deviance: 89.677  on 33  degrees of freedom
# Residual deviance: 34.307  on 31  degrees of freedom
# AIC: 1219
# Number of Fisher Scoring iterations: 1
# Theta:  18.48 
# Std. Err.:  4.44 
# 2 x log-likelihood:  -1211.019 

Anova(FitNB2)
# Analysis of Deviance Table (Type II tests)
# Response: Day10Spores
# LR Chisq Df Pr(>Chisq)    
# Treatment   55.371  2  9.472e-13 ***

summary(glht(FitNB2,linfct=mcp(Treatment="Tukey")))
# Linear Hypotheses:
#              Estimate Std. Error z value Pr(>|z|)    
# NP - B ==  -0.72006    0.09711  -7.415  < 1e-04 ***
# S - B ==   -0.39692    0.09920  -4.001 0.000182 ***
# S - NP ==   0.32314    0.09711   3.327 0.002524 ** 

anova(FitNB, FitNB2, test="Chisq")
# Likelihood ratio tests of Negative Binomial Models
# Response: Day15Spores
# Response: Day10Spores
#     Model   theta Resid.     df      2 x log-lik.         Test     df   LR stat.    Pr(Chi)
# 1 Treatment 32.75756         31       -1189.305                               
# 2 Treatment 18.47607         31       -1211.019          1 vs 2     0    -21.71399       1

Lsmeans.df2<-as.data.frame(summary(lsmeans(FitNB2, ~Treatment)))
Lsmeans.df2$btmean<-exp(Lsmeans.df2$lsmean)
Lsmeans.df2$btlo<-exp(Lsmeans.df2$lsmean - Lsmeans.df2$SE) #lower bound
Lsmeans.df2$bthi<-exp(Lsmeans.df2$lsmean + Lsmeans.df2$SE) #upper bound
View(Lsmeans.df2)
Lsmeans.df2
# Treatment   lsmean         SE df asymp.LCL asymp.UCL   btmean     btlo     bthi
# 1         B 18.24925 0.07014536 NA  18.11176  18.38673 84245455 78538524 90367074
# 2        NP 17.52918 0.06715907 NA  17.39755  17.66081 41004167 38340800 43852545
# 3         S 17.85232 0.07014536 NA  17.71484  17.98980 56645455 52808195 60761545

pollen.axis<-c( "Buckwheat","No-Pollen","Sunflower")
plot02<- ggplot(Lsmeans.df2, aes(x=Treatment, y=btmean)) + 
  geom_bar(position=position_dodge(), stat="identity", colour="black" ) +
  geom_errorbar(aes(ymin=btlo, ymax=bthi),
                size=1, width=0.2,  # Width of the error bars
                position=position_dodge(0.9))
plot02
plot012<-plot02 + theme_bw() #remove colored background
ylabel<- expression(bold(paste("Nosema spores/mL", sep="") )) 
plot022<-plot012 + ylab(ylabel) #y label
plot022
#relabel x axis:
plot052<-plot022 + xlab("Treatment") + # x axis label
  scale_x_discrete(breaks=levels(Lsmeans.df2$Treatment) , labels=pollen.axis)
plot052

plot062<-plot052 + theme(text=element_text (size=25, face="bold") )
plot062
##remove gridlines (option), remove border
neater2<- plot062 + theme(axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), 
                        panel.border = element_blank(),
                        panel.background = element_blank()) + theme(axis.line = element_line(size = 3, colour = "black"))
neater2


