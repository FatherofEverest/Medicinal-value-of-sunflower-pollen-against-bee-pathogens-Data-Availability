#### Sunflower Farm Descriptive stats with relation to Sun/Farm/Crop Area

rm(list=ls()) #clear memory

setwd("/Volumes/GoogleDrive/My Drive/POLLINATOR RESEARCH/PHD/Project_Sunflower/Analysis.Individual.Bee.and.Farm/Analyses/Farm Sampling")
farmdata<- read.csv("Analysis_1_All.csv", header = TRUE)
View(farmdata) 
str(farmdata)


##### Farm area ####

library(doBy)
FarmSizes<-summaryBy( FarmSize ~ Farm, data = farmdata, 
          FUN = function(x) { c(m = mean(x), s = sd(x)) } )

write.csv(FarmSizes, file = "Farm Sizes.csv")

farmdata2<-read.csv("Farm Sizes.csv", header = TRUE)
str(farmdata2)
mean(farmdata2$Area_m2) #163733.4 meters^2 = 0.163km^2


sd<-sd(farmdata2$Area_m2) #186781.6 m^2
SE<-sd/sqrt(22)
SE #39118.32 m^2 = 0.039km^2

range(farmdata2$Area_m2) #3458.35 m^2 to 628913.63 m^2 = 0.0035km^2 to 0.629km^2






mean(PropSunPerFarm) #0.009627425
sd(PropSunPerFarm) #0.02714787
range(PropSunPerFarm) #0.0000000 to 0.1387945




#### Crop area ####
PropSunPerCrop<- farmdata$SunArea/farmdata$Total.Crop.Area
PropSunPerCrop
farmdata$PropSunPerCrop<-PropSunPerCrop


mean(farmdata$Total.Crop.Area) #64639.57
sd(farmdata$Total.Crop.Area) #78268.6
range(farmdata$Total.Crop.Area) #675.0  to 295218.2

mean(PropSunPerCrop) #0.03435834
sd(PropSunPerCrop) #0.0519048
range(PropSunPerCrop) #0.0000000  to 0.2426264

PropCropPerFarm<-farmdata$Total.Crop.Area/farmdata$FarmSize
PropCropPerFarm
farmdata$PropCropPerFarm<-PropCropPerFarm

mean(PropCropPerFarm) #0.3285605
sd(PropCropPerFarm) #0.2425732
range(PropCropPerFarm) #0.005799247  to 1.000000000


cor(farmdata$SunArea, farmdata$FarmSize)
