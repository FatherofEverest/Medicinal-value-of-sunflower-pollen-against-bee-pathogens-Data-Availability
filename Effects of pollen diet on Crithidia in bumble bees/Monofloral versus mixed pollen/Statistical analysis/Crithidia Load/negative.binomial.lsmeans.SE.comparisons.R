###### Test of differences between glmmADMB, glmmTMB, and glmer.nb ##########
### Observed that glmmADMB confidence intervals were up to 8x larger than with other packages
## Test if this is consistent across: 
  ##1. Zero-inflated model
  ##2. Non-zero-inflated model
  ##3. Without random effects
  ##4. Poisson model

#I. Preliminaries:
#Load packages, load data, source function to extract lsmeans
detach
#install.packages("glmmADMB", repos = "http://glmmadmb.r-forge.r-project.org/repos")
#detach("package:Matrix", unload = TRUE)
library(glmmADMB)
library(glmmTMB)
library(lsmeans)
detach("package:Matrix", unload = TRUE)
#Use data from worked example
#http://glmmadmb.r-forge.r-project.org/glmmADMB.html

data(Owls)
str(Owls)
Owls <- transform(Owls, 
                  Nest=reorder(Nest,NegPerChick), 
                  logBroodSize=log(BroodSize), 
                  NCalls=SiblingNegotiation)

########   Interlude   #######
#Use Ben Bolker's function to talk to lsmeans
# https://github.com/glmmTMB/glmmTMB/issues/205
recover.data.glmmTMB <- function(object, ...) {
  fcall <- getCall(object)
  recover.data(fcall,delete.response(terms(object)),
               attr(model.frame(object),"na.action"), ...)
}
lsm.basis.glmmTMB <- function (object, trms, xlev, grid, vcov.,
                               mode = "asymptotic", component="cond", ...) {
  if (mode != "asymptotic") stop("only asymptotic mode is available")
  if (component != "cond") stop("only tested for conditional component")
  if (missing(vcov.)) 
    V <- as.matrix(vcov(object)[[component]])
  else V <- as.matrix(.my.vcov(object, vcov.))
  dfargs = misc = list()
  if (mode == "asymptotic") {
    dffun = function(k, dfargs) NA
  }
  ## use this? misc = .std.link.labels(family(object), misc)
  contrasts = attr(model.matrix(object), "contrasts")
  m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
  X = model.matrix(trms, m, contrasts.arg = contrasts)
  bhat = fixef(object)[[component]]
  if (length(bhat) < ncol(X)) {
    kept = match(names(bhat), dimnames(X)[[2]])
    bhat = NA * X[1, ]
    bhat[kept] = fixef(object)[[component]]
    modmat = model.matrix(trms, model.frame(object), contrasts.arg = contrasts)
    nbasis = estimability::nonest.basis(modmat)
  }
  else nbasis = estimability::all.estble
  list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
       dfargs = dfargs, misc = misc)
}

#####   End interlude ###



#######################################################
###1. With zero inflation ###########################
###################################################
m.zi<- glmmadmb(NCalls~FoodTreatment+ArrivalTime+ 
                  +(1|Nest), 
                data=Owls, 
                zeroInflation=TRUE, 
                family="nbinom")
summary(m.zi)
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)             4.2674     0.4705    9.07  < 2e-16 ***
#   FoodTreatmentSatiated  -0.2602     0.0845   -3.08   0.0021 ** 
#   ArrivalTime            -0.0840     0.0190   -4.42  9.8e-06 ***

#Plot lsmeans by FoodTreatment 
owls.lsm.zi<-lsmeans(m.zi, ~FoodTreatment)
owls.lsm.zi 
# FoodTreatment   lsmean        SE df asymp.LCL asymp.UCL
# Deprived      2.188727 0.7205142 NA 0.7765454  3.600909
# Satiated      1.928499 0.7498151 NA 0.4588887  3.398110

plot(owls.lsm.zi)
#95% confidence bands overlap almost entirely



##################  Compare glmmADMB fit to glmmTDMB  ####################
#install.packages("glmmTMB")
zi.t<- glmmTMB(NCalls~FoodTreatment+ArrivalTime+ 
                  +(1|Nest), 
                data=Owls, 
                ziformula = ~1,
                family="nbinom2")
summary(zi.t)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            4.26735    0.47044   9.071  < 2e-16 ***
#   FoodTreatmentSatiated -0.26022    0.08450  -3.080  0.00207 ** 
#   ArrivalTime           -0.08396    0.01898  -4.423 9.74e-06 ***

#Compare to glmmADMB model:Fixed effects are identical
summary(m.zi)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)             4.2674     0.4705    9.07  < 2e-16 ***
#   FoodTreatmentSatiated  -0.2602     0.0845   -3.08   0.0021 ** 
#   ArrivalTime            -0.0840     0.0190   -4.42  9.8e-06 ***

#Plot lsmeans by FoodTreatment 
####nb: Extract lsmeans from glmmTMB with the helper function at start of script ###

lsm.TMB<- lsmeans(zi.t, ~FoodTreatment)
plot(lsm.TMB)  #non-overlapping CI's

#Compare SE's
owls.lsm.zi #from ADMB
# FoodTreatment   lsmean        SE df asymp.LCL asymp.UCL
# Deprived      2.188727 0.7205142 NA 0.7765454  3.600909
# Satiated      1.928499 0.7498151 NA 0.4588887  3.398110


lsm.TMB
# FoodTreatment   lsmean         SE df asymp.LCL asymp.UCL
# Deprived      2.188720 0.06118962 NA  2.068790  2.308649
# Satiated      1.928498 0.08419132 NA  1.763486  2.093510

#lsmeans are identical but SE's differ by factor of 9 to 12?!


#######################################################
###2. Without zero inflation ###########################
###################################################

m.nb<- glmmadmb(NCalls~FoodTreatment+ArrivalTime+ 
                  +(1|Nest), 
                data=Owls, 
                zeroInflation=FALSE, 
                family="nbinom")
summary(m.nb)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)             4.9101     0.6334    7.75  9.1e-15 ***
#   FoodTreatmentSatiated  -0.6924     0.1069   -6.48  9.4e-11 ***
#   ArrivalTime            -0.1154     0.0253   -4.57  4.9e-06 ***
#Plot lsmeans by FoodTreatment 
owls.lsm<-lsmeans(m.nb, ~FoodTreatment)
owls.lsm 
# FoodTreatment   lsmean        SE df  asymp.LCL asymp.UCL
# Deprived      2.053073 0.8952071 NA  0.2984988  3.807646
# Satiated      1.360690 0.9037320 NA -0.4105918  3.131973

plot(owls.lsm)
#95% confidence bands overlap almost entirely

#Confirm with predict.glmmadmb:
New.data<-expand.grid(FoodTreatment= levels(Owls$FoodTreatment),
                      ArrivalTime = mean(Owls$ArrivalTime))

#Get standard errors:
calls.pred<- predict(m.nb, New.data, re.form = NA, se.fit = TRUE)
calls.pred<-data.frame(calls.pred)
New.data$NCalls <- calls.pred$fit
New.data$SE<-calls.pred$se.fit
New.data
# FoodTreatment ArrivalTime   NCalls        SE
# 1      Deprived    24.75763 2.053073 0.8952071
# 2      Satiated    24.75763 1.360690 0.9037320
#Matches with lsmeans output



##################  Compare glmmADMB fit to glmmTDMB  ####################
m.nb2<- glmmTMB(NCalls~FoodTreatment+ArrivalTime+ 
                  +(1|Nest), 
                data=Owls, 
                family="nbinom2")
summary(m.nb2)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            4.91011    0.63343   7.752 9.07e-15 ***
#   FoodTreatmentSatiated -0.69238    0.10692  -6.476 9.44e-11 ***
#   ArrivalTime           -0.11540    0.02526  -4.569 4.90e-06 ***

#Compare to glmmADMB model:Fixed effects are identical
summary(m.nb)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)             4.9101     0.6334    7.75  9.1e-15 ***
#   FoodTreatmentSatiated  -0.6924     0.1069   -6.48  9.4e-11 ***
#   ArrivalTime            -0.1154     0.0253   -4.57  4.9e-06 ***

#Plot lsmeans by FoodTreatment 
lsm.TMB<- lsmeans(m.nb2, ~FoodTreatment)
plot(lsm.TMB)  #non-overlapping CI's

#Compare SE's
owls.lsm
# FoodTreatment   lsmean        SE df  asymp.LCL asymp.UCL
# Deprived      2.053073 0.8952071 NA  0.2984988  3.807646
# Satiated      1.360690 0.9037320 NA -0.4105918  3.131973

lsm.TMB
# FoodTreatment   lsmean        SE df asymp.LCL asymp.UCL
# Deprived      2.053065 0.1068562 NA  1.843631  2.262500
# Satiated      1.360683 0.1161322 NA  1.133068  1.588298

#lsmeans are identical but SE's differ by factor of 8?!


##Compare to lme4 fit:
library(lme4)

nb4<-glmer.nb(NCalls~FoodTreatment+ArrivalTime+ 
                +(1|Nest), 
              data=Owls)

#Convergence warning
nb4.lsm<-lsmeans(nb4, ~FoodTreatment)

plot(nb4.lsm) #well-separated, glmmADMB SE's seem to be anomalous?

### Seems that differences in SE's are similar regardless of zero inflation


#######################################################
###3. Without random effects ###########################
###################################################

###ADMB
No.re.admb<- glmmadmb(NCalls~FoodTreatment+ArrivalTime, 
                      data=Owls, 
                      zeroInflation=FALSE, 
                      family="nbinom")

###TMB
No.re.tmb<- glmmTMB(NCalls~FoodTreatment+ArrivalTime, 
                    data=Owls, 
                    family="nbinom2")
  #Models fit quickly with no random effect !

#Confirm that fixed effects are identical:
fixef(No.re.admb)
# (Intercept) FoodTreatmentSatiated           ArrivalTime 
# 5.1002637            -0.4855420            -0.1221893 
fixef(No.re.tmb)
# (Intercept)  FoodTreatmentSatiated            ArrivalTime  
# 5.1002                -0.4855                -0.1222  

##lsmeans
lsm.a<-lsmeans(No.re.admb, ~FoodTreatment)
lsm.a
# FoodTreatment   lsmean        SE df  asymp.LCL asymp.UCL
# Deprived      2.075145 0.8920297 NA  0.3267991  3.823491
# Satiated      1.589603 0.8984394 NA -0.1713057  3.350512
lsm.t<-lsmeans(No.re.tmb, ~FoodTreatment)
lsm.t
# FoodTreatment   lsmean         SE df asymp.LCL asymp.UCL
# Deprived      2.075137 0.06685284 NA  1.944108  2.206166
# Satiated      1.589594 0.07356489 NA  1.445410  1.733779

#Here there is >10x difference in SE's !!

par(mfrow=c(2,1))
plot(lsm.a)
plot(lsm.t)
#quite different


## The irony here is that pairwise comparisons,
  #i.e., uncertainty for the DIFFERENCES... are the same ##
lsmeans(No.re.admb, pairwise~FoodTreatment)
# $contrasts
# contrast            estimate       SE df z.ratio p.value
# Deprived - Satiated 0.485542 0.099399 NA   4.885  <.0001

lsmeans(No.re.tmb, pairwise~FoodTreatment)
# $contrasts
# contrast             estimate         SE df z.ratio p.value
# Deprived - Satiated 0.4855427 0.09939888 NA   4.885  <.0001

#######################################################
###4. With poisson family model ###########################
###################################################

poi.a <- glmmadmb(NCalls~FoodTreatment+ArrivalTime + (1|Nest), 
                 data=Owls, 
                 zeroInflation=FALSE, 
                 family="poisson")

poi.t<- glmmTMB(NCalls ~FoodTreatment+ArrivalTime + (1|Nest), 
                    data=Owls, 
                family="poisson")

##Compare fixed effects and SE's
summary(poi.a)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            5.155380   0.245960   20.96   <2e-16 ***
#   FoodTreatmentSatiated -0.590390   0.035959  -16.42   <2e-16 ***
#   ArrivalTime           -0.129272   0.009261  -13.96   <2e-16 ***
summary(poi.t)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            5.32204    0.25689    20.7   <2e-16 ***
#   FoodTreatmentSatiated -0.66440    0.03726   -17.8   <2e-16 ***
#   ArrivalTime           -0.13642    0.00962   -14.2   <2e-16 ***

#Estimates and SE's slightly different but not vastly different


(lsm.poi.a<- lsmeans(poi.a, pairwise ~ FoodTreatment))
# FoodTreatment   lsmean        SE df asymp.LCL asymp.UCL
# Deprived      1.944592 0.3577493 NA 1.2434164  2.645768
# Satiated      1.280187 0.3607572 NA 0.5731162  1.987258

# contrast             estimate      SE df z.ratio p.value
# Deprived - Satiated 0.6644048 0.03726 NA  17.832  <.0001

(lsm.poi.t<- lsmeans(poi.t, pairwise ~ FoodTreatment))
# FoodTreatment   lsmean         SE df asymp.LCL asymp.UCL
# Deprived      1.954918 0.09723980 NA  1.764332  2.145505
# Satiated      1.364529 0.09908944 NA  1.170317  1.558740

# contrast             estimate        SE df z.ratio p.value
# Deprived - Satiated 0.5903898 0.0359587 NA  16.419  <.0001

plot(lsm.poi.a$lsmeans)
plot(lsm.poi.t$lsmeans)

###Same pattern as for negative binomial. 
##lsmeans SE's are 3x bigger for ADMB, but contrast SE's are similar


#Compare to glmer:
poi.lme <- glmer(NCalls~FoodTreatment+ArrivalTime + (1|Nest), 
                  data=Owls, 
                  family="poisson")
#plot all three:
par(mfrow=c(1,3))
plot(lsm.poi.a$lsmeans)
plot(lsm.poi.t$lsmeans)
plot(lsmeans(poi.lme, pairwise~FoodTreatment)$lsmeans)
#Again the glmmadmb estimates differ from those for glmmTMB and glmer


#try mcmcglmm
library(MCMCglmm)

MCp<-MCMCglmm(fixed = NCalls~FoodTreatment+ArrivalTime,
         random = ~Nest,
         data = Owls,
         family = "poisson")
summary(MCp)
lsm.mcp<-lsmeans(MCp, ~FoodTreatment, data = Owls)
lsm.mcp
# FoodTreatment    lsmean        SE df asymp.LCL asymp.UCL
# Deprived      1.6753579 0.1494447 NA 1.3824516 1.9682641
# Satiated      0.6345644 0.1576183 NA 0.3256381 0.9434906
plot(lsm.mcp)

lsm.poi.t
# FoodTreatment   lsmean         SE df asymp.LCL asymp.UCL
# Deprived      1.954918 0.09723980 NA  1.764332  2.145505
# Satiated      1.364529 0.09908944 NA  1.170317  1.558740

poi.lme

lsmeans(poi.lme, ~FoodTreatment)
# FoodTreatment   lsmean         SE df asymp.LCL asymp.UCL
# Deprived      1.954931 0.09719423 NA  1.764434  2.145428
# Satiated      1.364542 0.09903617 NA  1.170435  1.558650



plot(lsmeans(poi.lme, ~FoodTreatment))
#similar
plot(lsmeans(poi.a, ~FoodTreatment))
#way bigger CI's
plot(lsmeans(poi.t, ~FoodTreatment))


#SE's are in same ballpark for both MCMC, TMB, and lme4, but bigger from glmmADMB
#MCMC predicted  much lower mean for Satiated 
#compare to actual mean
library(plyr)
library(Rmisc)
owls.sum<- summarySE(data = Owls, measurevar = "NCalls", groupvars = "FoodTreatment" )
owls.sum$logCalls<- log(owls.sum$NCalls)
owls.sum

# FoodTreatment   N   NCalls       sd        se        ci logCalls
# 1      Deprived 320 8.162500 6.454458 0.3608152 0.7098780 2.099550
# 2      Satiated 279 5.064516 6.540092 0.3915451 0.7707697 1.622259
 #MCMCglmm estimates are quite low (mu = 0.63) for Satiated 
 #This round goes to the frequentist models!

