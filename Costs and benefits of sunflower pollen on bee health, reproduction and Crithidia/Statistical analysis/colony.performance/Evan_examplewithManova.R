ls()
rm(list=ls())
ls()
#Import Data set in R
alkaloid <- read.csv("~/Desktop/Manuscript /Chapter 4/Stats and data/AlkaloidParasite.csv")


#View Dataset and assign factors 
str(alkaloid)
"alkaloid$Height <- as.numeric(alkaloid$Height)
alkaloid$Mass <- as.numeric(alkaloid$Mass)"



#Load 'nlme' package
library("nlme", lib.loc=
          "/Library/Frameworks/R.framework/Versions/2.12/Resources/library")

source('~//Dropbox/muvari/biostats.R')


##checking normality

qqnorm(alkaloid$Xlupanine)
shapiro.test(alkaloid$Xlupanine)

qqnorm(alkaloid$Lupanine)
shapiro.test(alkaloid$Lupanine)
qqnorm(alkaloid$Oxolupanine)
shapiro.test(alkaloid$Oxolupanine)
qqnorm(alkaloid$Hydroxylupanine)
shapiro.test(alkaloid$Hydroxylupanine)

##############

"alkaloid2<-alkaloid[c(-5,-7,-10,-8),]"

###Running Manova
alkaloid2<-cbind(alkaloid$Lupanine,alkaloid$Xlupanine,alkaloid$Oxolupanine,alkaloid$Hydroxylupanine)

Manova<-manova(alkaloid2~alkaloid$Treat)
Manova
summary(Manova,test="Pillai")


###Running model and anova for alkaloids

model1<-lm(Lupanine~Treat,data=alkaloid)
model1
anova(model1)
summary(fml<-aov(Lupanine~Treat, data=alkaloid))


#######################################

model2<-lm(Xlupanine~Treat, data=alkaloid)
model2
anova(model2)
summary(fml<-aov(Xlupanine~Treat, data=alkaloid))

##################################
model3<-lm(Hydroxylupanine~Treat, data=alkaloid)
model3
anova(model3)
summary(fml<-aov(Hydroxylupanine~Treat, data=alkaloid))

#########################
model4<-lm(Oxolupanine~Treat,data=alkaloid)
model4
anova(model4)
summary(fml<-aov(Oxolupanine~Treat, data=alkaloid))

###################means and SEs

summarySE(alkaloid, measurevar="Lupanine", groupvars=c("Treat"))
summarySE(alkaloid, measurevar="Xlupanine", groupvars=c("Treat"))
summarySE(alkaloid, measurevar="Oxolupanine", groupvars=c("Treat"))
summarySE(alkaloid, measurevar="Hydroxylupanine", groupvars=c("Treat"))

##################checking for residuals
no.tf.anova<-lm(Xlupanine~Treat, data=alkaloid2)
summary(no.tf.anova)
plot(no.tf.anova)
plot(no.tf.anova$residuals)

no.tf.anova<-lm(Hydroxylupanine~Treat, data=alkaloid2)
summary(no.tf.anova)
plot(no.tf.anova)
plot(no.tf.anova$residuals)

no.tf.anova<-lm(Lupanine~Treat, data=alkaloid2)
summary(no.tf.anova)
plot(no.tf.anova)
plot(no.tf.anova$residuals)



