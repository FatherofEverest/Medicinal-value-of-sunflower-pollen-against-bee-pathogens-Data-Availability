
#install.packages("glmmADMB", repos = "http://glmmadmb.r-forge.r-project.org/repos")

library(glmmADMB)
library(glmmTMB)
library(lsmeans)

#Use data from worked example
#http://glmmadmb.r-forge.r-project.org/glmmADMB.html

data(Owls)
str(Owls)
Owls <- transform(Owls, 
                  Nest=reorder(Nest,NegPerChick), 
                  logBroodSize=log(BroodSize), 
                  NCalls=SiblingNegotiation)


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
#SE is much higher than for fixed effects in model

plot(owls.lsm)
 #95% confidence bands overlap almost entirely

#Confirm with predict.glmmadmb:
New.data<-expand.grid(FoodTreatment= levels(Owls$FoodTreatment),
                      ArrivalTime = mean(Owls$ArrivalTime))

New.data$NCalls <- predict(m.nb, New.data, re.form=NA, SE.fit = TRUE)

#Get standard errors:
calls.pred<- predict(m.nb, New.data, re.form = NA, se.fit = TRUE)
calls.pred<-data.frame(calls.pred)

New.data$SE<-calls.pred$se.fit
New.data
# FoodTreatment ArrivalTime   NCalls        SE
# 1      Deprived    24.75763 2.053073 0.8952071
# 2      Satiated    24.75763 1.360690 0.9037320
#Matches with lsmeans output



##################  Compare to glmmTDMB  ####################
#install.packages("glmmTMB")
library(glmmTMB)
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
owls.lsm<-lsmeans(m.nb2, ~FoodTreatment)
#oops, lsmeans can't use glmmTMB object!

  ########   Interlude   #######
#Ben Bolker wrote a function to talk to lsmeans-- incredible!
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

#Confirm with predict()
New.data2<-expand.grid(FoodTreatment= levels(Owls$FoodTreatment),
                      ArrivalTime = mean(Owls$ArrivalTime) ) #,Nest = NA)

New.data2$NCalls <- predict(m.nb2, New.data2, re.form=NA, SE.fit = TRUE)
# Error in predict.glmmTMB(m.nb2, New.data2, re.form = NA, SE.fit = TRUE) : 
#   re.form not yet implemented
predict(m.nb2, newdata = New.data2)
#Error in eval(predvars, data, env) : object 'Nest' not found
New.data2$Nest <- "Forel"
predict(m.nb2, New.data2, se.fit = TRUE)


##Reply from John Maindonald:
 #You are getting estimates for contrasts, not for the groups
# The confidence intervals that you have obtained are for the levels
# of `FoodTreatment`, not for the contrast `Satiated-Deprived`.  
 
# Try, the following, which also gives a confidence interval for the
# difference from the initial level of `FoodTreatment`:
  
library(glmmADMB)
library(lsmeans)
 Owls <- transform(Owls,
                    Nest=reorder(Nest,NegPerChick),
                    logBroodSize=log(BroodSize),
                    NCalls=SiblingNegotiation)
 m.nb<- glmmadmb(NCalls~FoodTreatment + ArrivalTime + (1|Nest),
                 data = Owls, zeroInflation=TRUE, family = "nbinom")
   
   owls.lsm<-lsmeans(m.nb, ~FoodTreatment)
    lsmeans (owls.lsm, "FoodTreatment", contr = "trt.vs.ctrl")
                  # $lsmeans
                  # . . .
                  # 
                  # $contrasts
                  # contrast             estimate       SE df z.ratio p.value
                  # Satiated - Deprived -0.260228 0.084501 NA   -3.08  0.0021
    plot(owls.lsm)

#still overlapping
    
    
##Second reply:
    # The differences between the two graphs are a conseqence of
    # what is done on the way to creating those graphs.  For comparing
    # levels of `FoodTreatment`, the SE is for the
    # difference, not for the levels individually.  It is then information that it
    # would be helpful to add to the graph given by
  
    owls.lsm<-lsmeans(m.nb, ~FoodTreatment)
    plot(owls.lsm)
    
## One can get a plot that shows the SE for the difference thus:
    coef(m.nb)
K <- diag(length(coef(m.nb)))[1:2,]
K    #this will contrast foodtreatment Deprived vs foodtreatment satiated
rownames(K) <- c("Deprived", "Sat-Dep")
library(multcomp)
glht(m.nb, linfct = K)
# Estimate
# Deprived == 0   4.9101
# Sat-Dep == 0   -0.6924
plot(glht(m.nb,linfct=K))
                     
# Or, nearer to what you want, maybe:
K2 <- rbind(K[1,], c(1,1,0), K[2,])
K2
#Deprived; Satiated; Satiated - Deprived?
#Second row looks strange, how can we have Satiated and deprived together?
rownames(K2) <- c("Deprived","Saturated", "Sat-Dep")
plot(glht(m.nb,linfct=K2)) #??what is the second row? 
  #How can nest be both satiated and deprived? 
                     # It is, of course, in this simple case, possible to place intervals
                     # around the two estimates, designed so that if the intervals do
                     # not overlap, then the difference is not "significant" at alpha=0.05.
                     # 
                     # I will leave it to you, or to others, to check just what the code you
                     # give, that uses as its starting-point output from glmmTMB(), may be
                     # doing.  This is not a straightforward use of lsmeans().
                     # 
                     # John Maindonald             email: john.maindonald@anu.edu.a 
                     # 
  
    
##Compare to lme4 fit:
library(lme4)
nb4<-glmer.nb(NCalls~FoodTreatment+ArrivalTime+ 
                +(1|Nest), 
              data=Owls)
#Convergence warning
nb4.lsm<-lsmeans(nb4, ~FoodTreatment)
plot(nb4.lsm) #well-separated


##Check Russell Lenth's function to extract data from glmmADMB
#http://permalink.gmane.org/gmane.comp.lang.r.lme4.devel/12053
#Similar to Ben Bolker's function to extract from glmmTMB
#https://github.com/glmmTMB/glmmTMB/issues/205
recover.data.glmmadmb = lsmeans:::recover.data.lm

lsm.basis.glmmadmb = function (object, trms, xlev, grid)
{
  contrasts = object$contrasts
  m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
  X = model.matrix(trms, m, contrasts.arg = contrasts)
  bhat = fixef(object)
  V = vcov(object)
  misc = list()
  if (!is.null(object$family)) {
    fam = object$family
    misc$tran = object$link
    misc$inv.lbl = "response"
    if (!is.na(pmatch(fam,"binomial")))
      misc$inv.lbl = "prob"
    else if (!is.na(pmatch(fam,"poisson")))
      misc$inv.lbl = "rate"
  }
  nbasis = matrix(NA)
  dffun = function(...) NA
  list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun,
       dfargs = list(), misc = misc)
}
plot(lsmeans(m.nb, ~FoodTreatment))
#no change