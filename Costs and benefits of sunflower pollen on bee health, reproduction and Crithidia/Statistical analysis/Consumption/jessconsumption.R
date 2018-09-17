#Analysis of Jess Leslie nectar/pollen limitation experiment
#Measuring nectar consumption

#1. Prelim analysis on raw data (not accounting for evaporation)
#2. Second analysis will account for evaporation on each day

rm(list=ls())
#set wd
setwd("C:/Users/lynnsa/Documents/Dropbox/Adler Research Projects/Bee performance/Jess Sum15/jess stats")
#read in consumption data
Consumption<-read.csv('ConsumptionDataR.csv')
#read in data frame that has mass & colony for each bee
BeeData<-read.csv('C:/Users/Evan/Google Drive/bee/conroy.taylor/Data & Results/survival/Survival Data V2.csv')
#column for number of days until inoculation
setwd("C:/Users/Evan/Google Drive/bee/conroy.taylor/Data & Results")
Callows<-read.csv('CallowDataR.csv', na.strings = "NA")
#add the number of days from emergence to inoculation
BeeData$Num_of_Days= Callows$Num_of_Days
#change back wd
setwd("C:/Users/Evan/Google Drive/bee/conroy.taylor/Data & Results/consumption")
#merge BeeData with Consumption
#View(BeeData)
#View(Consumption)
library(plyr)
Cons.plus<-join(x=BeeData, y=Consumption, by="ID")
#View(Cons.plus) #great, new data frame with full data for each bee

#now integrate data from Evaporation controls
#second data frame for evaporation
#Has average values from the 6 reps of each sugar conc on each day
#column 30_avg_evap for 30% sucrose;
#15_avg_evap for 15% sucrose
#impute mean evaporation value for days where average evap is missing
controls<-read.csv('EvaporationControlsR.csv')
#View(controls)
controls.lowsugar<-subset(controls, !X15_Avg_Evap=="NA")
controls.lowsugar$X15_Avg_Evap
lowsugarmean<-mean(controls.lowsugar$X15_Avg_Evap) #118 mg
controls.highsugar<-subset(controls, !X30_Avg_Evap=="NA")
controls.highsugar$X30_Avg_Evap
highsugarmean<-mean(controls.highsugar$X30_Avg_Evap) 
highsugarmean ##108 mg
###this is a considerable difference!!!
#I feel that using these controls is going to throw the analysis
#as compared to using the raw consumption data from the bees' vials
#this feels really dangerous, I do this more as a coding exercise
Controls.imputed<-as.data.frame(controls)
#View(Controls.imputed)
(Controls.imputed$X15_Avg_Evap) 
df=Controls.imputed
str(df$X15_Avg_Evap)
mean(df$X15_Avg_Evap, na.rm=T)
df$X15_Avg_Evap[is.na(df$X15_Avg_Evap) == "TRUE"] <- mean(df$X15_Avg_Evap, na.rm = T)
df$X15_Avg_Evap #substituted means
df$X30_Avg_Evap[is.na(df$X30_Avg_Evap) == "TRUE"] <- mean(df$X30_Avg_Evap, na.rm = T)
df$X30_Avg_Evap #means are substituted
Controls.imputed=df
View(Controls.imputed)
#need to unstack data ? how to merge with consumption data?
#want 2 rows for each date: 
#one with 30% treat's evap avg; one with 15% treat's evap avg
library(reshape2)
#convert to long format
#http://www.cookbook-r.com/Manipulating_data/Converting_data_between_wide_and_long_format/
Imputed_long <- melt(Controls.imputed, 
                  # ID variables - all the variables to keep but not split apart on
                  id.vars=c("Date"),
                  # The source columns
                  measure.vars=c("X30_Avg_Evap", "X15_Avg_Evap" ),
                  # Name of the destination column that will identify the original
                  # column that the measurement came from
                  variable.name="Nec_Treat",
                  value.name="Avg_Evap"
)
#View(Imputed_long) 
#Average evaporation on each date for each concentration
# Rename factor names to match levels of Nec_Treat in bee data ("Cons.plus")
levels(Imputed_long$Nec_Treat)[levels(Imputed_long$Nec_Treat)=="X30_Avg_Evap"] <- "High"
levels(Imputed_long$Nec_Treat)[levels(Imputed_long$Nec_Treat)=="X15_Avg_Evap"] <- "Low"
#View(Imputed_long) #looks good

conctest<-aov(Avg_Evap~Nec_Treat, data=Imputed_long)
summary(conctest) #p=0.3, non-significant effect of nectar treatment on evaporation

#now merge the data 
Cons.complete<-join(x=Cons.plus, y=Imputed_long,
                     by=c("Date", "Nec_Treat"))
#View(Cons.complete) #rats, we are still missing some dates

#Let's use the means from the evaporation controls again
#ex df$Items[which(df$Store.Type == "A" | df$Store.Type == "C" )] <- 0
df<-Cons.complete
df$Avg_Evap[which(df$Nec_Treat == "High")] <- highsugarmean
df$Avg_Evap[which(df$Nec_Treat == "Low")] <- lowsugarmean
#View(df)
Data<-df

#Ah, now we finally are ready to look at the data and run the mixed model!
library(lme4)
#check the data-- quality controls; uncomment to view
#View(Data)
#Omit anything >0.5 g
#because bees were only feed 500 uL nectar!
Toomuch<-subset(Data, N_Consumed_g>=0.5)
#View(Toomuch)

RealConsumption<-subset(Data, N_Consumed_g<=0.5)
hist(RealConsumption$N_Consumed_g) #perfect distribution there
#New column for net consumption
RealConsumption$Nectar_Net<-RealConsumption$N_Consumed_g - RealConsumption$Avg_Evap
#quick plot
plot(RealConsumption$Nectar_Net~RealConsumption$Group)
#media 180 mg before evaporation correction, now <100 mg
#looks like less consumption in nectar-limited groups?
#possibly an artifact of differences betweeen 
#highsugar and lowsugar controls evaporation values

RealConsumption$ID<-as.factor(RealConsumption$ID)
RealConsumption$Day_Num<-as.factor(RealConsumption$Day_Num)

testmodel<-aov(Nectar_Net~Pol_Treatment*Nec_Treatment, data=RealConsumption)
summary(testmodel) #with corrected value, nectar has an effect
#properly we should use mixed model to account for repeated measures:
library(lme4)
nectar.2rand0 <- lmer(
  Nectar_Net~Pol_Treatment*Nec_Treatment + Colony + 
    Callow_Mass_g + Num_of_Days + 
    (1|Day_Num) + (1|ID) + (1|Date) , data=RealConsumption) 
summary(nectar.2rand0)
library(car)
Anova(nectar.2rand0) #eliminate "Num_of_Days"

nectar.2rand <- lmer(
  Nectar_Net~Pol_Treatment*Nec_Treatment + Colony + Callow_Mass_g +
    (1|Day_Num) + (1|ID) + (1|Date), data=RealConsumption) 
summary(nectar.2rand)
library(car)
Anova(nectar.2rand) #eliminate callow mass
plot(nectar.2rand) #these seem OK

#omit callow mass
NoCallow <- lmer(
  Nectar_Net~Pol_Treatment*Nec_Treatment + Colony  +
    (1|Day_Num) + (1|ID) +(1|Date), data=RealConsumption) 
summary(NoCallow)
library(car)
Anova(NoCallow) 
plot(NoCallow) #Good

#try without days since inoculation 
#(remove the 'Day_Num' random effect)
NoCallow.NoTime<-lmer(
  Nectar_Net~Pol_Treatment*Nec_Treatment + Colony  +
     + (1|ID) + (1|Date), data=RealConsumption) 
summary(NoCallow.NoTime)
AIC(NoCallow) #-1206
AIC(NoCallow.NoTime) #-1151; not as good; keep "Day_Num"
anova(NoCallow, NoCallow.NoTime) #p=0.14, keep "Day_Num"
#we will retain the 'Day_Num' random effect

#test effect of experimental date:
NoCallownodate <- lmer(
  Nectar_Net~Pol_Treatment*Nec_Treatment + Colony  +
    (1|Day_Num) + (1|ID) , data=RealConsumption) 
summary(NoCallownodate)
anova(NoCallow, NoCallownodate)
#highly significant effect of date; retain in model
#Report "NoCallow" as reduced model
RedMod<-NoCallow
Anova(RedMod)
library(lsmeans)
lsmeans(NoCallow, ~Pol_Treatment*Nec_Treatment)
lsmeans(NoCallow, ~Nec_Treatment)
#about 15 mg higher consumption in High sugar treatments

#for comparison, here is the model without correction for evaporation
Nocontrols<-lmer(
  N_Consumed_g~Pol_Treatment*Nec_Treatment + Colony  +
    + (1|Day_Num) + (1|ID) + (1|Date) , data=RealConsumption) 
Anova(Nocontrols) #now there is no effect of the nectar treatment
lsmeans(Nocontrols, ~Pol_Treatment*Nec_Treatment)

#I would report "NoCallow" as the final model


My.Lsmeans<-(lsmeans(NoCallow, ~Pol_Treatment*Nec_Treatment))
My.Lsmeans
#######PLOTTING #############
#plotting marginal means:
library(ggplot2)
#make a data frame from the lsmeans object
Lsmeans.df<-summary(My.Lsmeans)
str(Lsmeans.df)
View(Lsmeans.df)
summary(Lsmeans.df$Nec_Treatment)
#rename sugar concentrations to "High" and "Low"
library(plyr)
Lsmeans.df$Nectar.conc<-revalue(Lsmeans.df$Nec_Treatment, 
        c("HighNectar"="High", "LowNectar"="Low"))
#rename column for compatibility with previous script
Lsmeans.df$Pollen<-Lsmeans.df$Pol_Treatment
View(Lsmeans.df)
#reorder the levels of pollen treatment
Lsmeans.df$Pollen=factor(Lsmeans.df$Pollen, 
               levels (Lsmeans.df$Pollen) [c(2,1)] )
levels(Lsmeans.df$Nectar.conc)
levels(Lsmeans.df$Pollen)
#label for x axis
pollen.axis<-c("Pollen", "No Pollen")
library(ggplot2)
plot0<- ggplot(Lsmeans.df, aes(x=Pol_Treatment, y=lsmean, fill=Nectar.conc)) + 
  geom_bar(position=position_dodge(), stat="identity", colour="black") +
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE),
                size=1, width=.2,                    # Width of the error bars
                position=position_dodge(.9))
plot0 #Solid base plot
plot01<-plot0 + theme_bw() #remove colored background
#superscripts for axis
#expression(paste('(SOC[',kgm^-2,'])'^0.25,sep=''))
ylabel<-expression(bold(paste ('Nectar consumption (g*', day^-1, ')')))
plot02<-plot01 + ylab(ylabel) #y label
#for slash instead of superscript:
#plot02<-plot01 + ylab("Nectar consumption (g/day)") #y label
plot03<- plot02 + scale_fill_hue(name="Sugar \nConcentration") +
  scale_x_discrete(breaks=levels(Lsmeans.df$Pollen) , labels=pollen.axis)
plot03 
plot04<- plot03 +
  scale_fill_manual(name="Sugar \nConcentration", #change legend name
                    breaks=c("High", "Low"), 
                    values=c("gray", "white"))#change fill
plot04
plot004<- plot04 + guides(fill = 
                            guide_legend(override.aes = list(colour = NULL))) +
  theme(legend.key=element_rect(colour = 'black', size=1))
plot004 #no more legend slashes, but kept outlines 

#relabel x axis:
plot05<-plot004 + xlab(NULL)  #removes x axis label
plot05
#Larger font:
plot06<-plot05 + theme(text=element_text (size=25, face="bold") )
plot06
##remove gridlines, remove border
neater<- plot06 +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.line = element_line(size = 1.5, colour = "black"))
neater
#justify axis to start at 0
Zeroed<- neater + coord_cartesian(ylim = c(0, 0.123))
Zeroed
#move legend
MoveLeg<-Zeroed + theme(legend.position=c(0.5, 0.9),
  legend.background = element_rect(color="black")) #box around legend
MoveLeg
#if we want y axis to start at 0:
#this removes space between x axis and bottom of bar
#Justified<- neater +
#scale_y_continuous(expand = c(0,0), limits = c(0,0.120))
#Justified
#adjust legend
fig3<-MoveLeg #or "Justified" to start y axis at 0

#export to PDF
pdf("fig3_consumption.pdf",height=6,width=6,paper='special') 
fig3
dev.off()

#export to eps
setEPS()
postscript("fig3nectarconsumption.eps", horizontal = FALSE, 
           onefile = FALSE, paper = "special", 
           height = 6, width = 6)
fig3
dev.off()

