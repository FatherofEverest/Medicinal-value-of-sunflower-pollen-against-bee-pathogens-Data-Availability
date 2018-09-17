######SCRIPT TO ANALYZE CRITHIDIA FROM WILD BEES########
###Collected from Amherst area farms summer 2015
####PROPOSED EXPLANATORY VARIABLE: sunflower area#

#libraries
#most of these we don't actually need here (epy)
library(MASS)
library(car)
library(glmmADMB)
library(bbmle)
library(effects)
library(lsmeans)
library(ggplot2)
#install.packages("geoR")
library(geoR)
#install.packages("ape")
library(ape)
#install.packages("ade4")
library(ade4)

#read & process data

farmdata<-read.csv("C:/Users/Evan/Google Drive/Project_Sunflower/Analysis.Individual.Bee.and.Farm/Analyses/Farm Sampling/Analysis_1_All.csv")
  #("C:/Users/csutherland/Dropbox/Adler/adler_bee.csv",header=TRUE)

coord <- read.table("C:/Users/Evan/Google Drive/Project_Sunflower/Analysis.Individual.Bee.and.Farm/Analyses/Farm Sampling/Sutherland_spatial_issues/coords.txt",
                    header=TRUE)
#("C:/Users/csutherland/Dropbox/Adler/coords.txt",header=TRUE)
str(coord)
length(farmdata$Farm) #667
#check that the levels match: (they don't)
levels(farmdata$Farm)
levels(coord$Farm)
#manually change the farm names in "coord" to match those in "farmdata"
library(plyr)
coord$Farm<-revalue(coord$Farm, c("LAU"="LAR", "DNG"="DAV"))

farmdataplus <- merge(farmdata,coord,by="Farm")
length(farmdataplus$Farm) #667, good

d<- data.frame(farmdata$CollectionDate,farmdata$Log1P.SunArea.,farmdata$WingSize)

farmdata$Count<-as.numeric(farmdata$Count)

#######RUN MODEL to generate residuals######
L1plus<-glmmadmb(Count ~ Log1P.SunArea. + CollectionDate +(1|Farm),
                 family = "nbinom",zeroInflation = FALSE, data=farmdata)
residuals(L1plus)
farmdataplus$X
#visualize residuals
sac <- cbind(farmdataplus$X,farmdataplus$Y,residuals(L1plus))
str(sac)
#View(sac) #finally got the X and Y columns, whew
set.seed(1)

par(mar=c(4,4,0,0),oma=c(0,0,0,0))


#plot map of residuals and export
setwd("C:/Users/Evan/Google Drive/Project_Sunflower/Analysis.Individual.Bee.and.Farm/Analyses/Farm Sampling/Sutherland_spatial_issues")
pdf("cs_variogram.pdf")
par(mar=c(5,5,4,2)+0.1,mgp=c(3,1,0))
plot(jitter(sac[,c(2,1)]),cex=sqrt(1+sac[,3]),pch=16,
     col=adjustcolor(ifelse(sac[,3]>0,4,2),alpha.f=0.3),
     xlab=expression(~Longitude~(degree)), 
     ylab=expression(~Latitude~(degree)))

dev.off()
#example of degrees in axis:
#ylab=expression(~degree~C)

#Moran's I
dmat <- 1/as.matrix(dist(sac[,c(2,1)]))

dmat[is.infinite(dmat)] <- 0

diag(dmat) <- 0

Moran.I(sac[,3], dmat)
#$observed
#[1] -0.004229671
#$expected
#[1] -0.001501502
#$sd
#[1] 0.00187081
#$p.value
#[1] 0.1447626
#FAIL to reject the null hypothesis that there is zero spatial
# autocorrelation present in the residuals :)


#Mantel test
dmat <- dist(sac[,c(2,1)])
rmat <- dist(sac[,3])
mantel.rtest(dmat, rmat, nrepet = 9999) #takes some time
#FAIL to reject the null hypothesis that there is zero spatial
# autocorrelation present in the residuals :)
#> mantel.rtest(dmat, rmat, nrepet = 9999)
#Monte-Carlo test
#Monte-Carlo test
#Observation: -0.01246183 
#Call: mantel.rtest(m1 = dmat, m2 = rmat, nrepet = 9999)
#Based on 9999 replicates
#Simulated p-value: 0.67 
#Warning messages:
#  1: In is.euclid(m1) : Zero distance(s)
#2: In is.euclid(m2) : Zero distance(s)
#3: In is.euclid(distmat) : Zero distance(s)

#(semi)variogram - should be flat and 
#high valued across distance bins
# meaning no correlation structure 
#THIS IS FLAT ACROSS THE RANGE
dmat <- dist(sac[,c(2,1)])
maxD <- max(dmat)
vg <- variog(coords = sac[,1:2], data = sac[,3],
             breaks = seq(0, maxD, l = 11))
#pdf("cs_semi_variogram.pdf") #uncomment to export
plot(vg, type = "b", main = "Variogram: residuals")
#dev.off()


