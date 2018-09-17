#script to convert sunflower microcolony data from long to wide format
#use "melt" function in "plyr"

#Goal:
#retain coluns for microcolony, infection, pollen
#one column for sum of larval masses by microcolony
#one column for whether or not pupae produced
#one column for number of eggs
#one column for number of honey pots
#one column for number of dead larvae

Df<-read.csv('C:/Users/Evan/Dropbox/Jess Sum15/Pollen consumption trials/performance.160228.csv',
            na.strings = ".")
View(Df)
#we need to have rows with counts for larval, pupal, callow
#right now we are missing them, with "stage-NA"
#For now, recode "stage=NA" to stage=pupal

Df<-Df [ 1:213, 1:8]
#View(Df)
mydata<-Df

mydata[,c("Mass", "Honeypots", "Count.dead")] <-
  apply(mydata[,c("Mass", "Honeypots", "Count.dead")], 2, function(x){replace(x, is.na(x), 0)})
View(mydata)
str(mydata$Stage)
mydata$Stage

mydata$Stage<-replace(mydata$Stage, is.na(mydata$Stage), "larval") 
mydata$Stage
View(mydata)

#make a new column for mass of just larvae
mydata$MassL<- mydata$Mass
mydata$MassL[mydata$Stage  != "larval"] <- 0
View(mydata)

#first we need to split apart based on "stage"
library(plyr)
library(reshape2)

Df<-mydata

wider<-dcast(Df, Microcolony + Infection + Pollen + Mass + Honeypots + Count.dead + MassL~ Stage, value.var="Count")
View(wider)

#substitute "0" for NA to avoid dropping microcolonies
wider[is.na(wider)] <- 0

#now summarize by microcolony
library(plyr)
#we need to sum across mass, honeypots, count.dead, callow, egg, larvae, pupae, NA (meaning callow and larvae and pupae?)

SUMMARIZED<-ddply(wider, c("Microcolony"), summarize,
                                TOTMASS=sum(Mass, na.rm=TRUE),
                                HPOTS=sum(Honeypots, na.rm=TRUE),
                                NUMDEAD=sum(Count.dead, na.rm=TRUE),
                                NUMCALLOWS=sum(callow, na.rm=TRUE),
                                NUMEGGS=sum(egg, na.rm=TRUE),
                                NUMLARV=sum(larval, na.rm=TRUE),
                                NUMPUP=sum(pupal, na.rm=TRUE),
                                LARVALMASS=sum(MassL))
View(SUMMARIZED)
#$TOTMASS : here I added up mass for eggs + larvael + pupae

#note-- no microcolonies 50 or 54
#now make a column that includes callows with pupae
SUMMARIZED$Totpup<- SUMMARIZED$NUMCALLOWS + SUMMARIZED$NUMPUP


#or merge back with plyr
summaryplus<-join(x=SUMMARIZED, y=wider, by="Microcolony", match="first")
View(summaryplus)
summaryplus2<- summaryplus[,1:11]
View(summaryplus2)
#missing micro's 50, 54-- these were removed beforehand
#also note that "TOTMASS is mass of ALL reproductive units (callow + pupae + egg)
#nice! 

#finally add back the source colonies
DF<-summaryplus2
DF$Colony<-DF$Microcolony
DF$Colony[DF$Microcolony  <=20] <- "A"
DF$Colony
DF$Colony[DF$Microcolony  >=21 & DF$Microcolony<=40 ] <- "B"
DF$Colony
DF$Colony[DF$Microcolony  >=41 & DF$Microcolony<=60 ] <- "C"
DF$Colony[DF$Microcolony  >=61 & DF$Microcolony<=80 ] <- "D"
DF$Colony

setwd("C:/Users/Evan/Dropbox/Jess Sum15/jess stats/colony.performance")
#write out a csv file for subsequent analysis! Nice job!
#Giving up is for the .... people who want to save time & be practical!
write.csv(DF, file="colony.performance.epy.practice.csv")
#remember to omit colonies 12 and 16 when analyzing!!
