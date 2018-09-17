
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



##################  Compare glmmADMB fit to glmmTDMB  ####################
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

#### Extract lsmeans from glmmTMB with the helper function ###

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
                      ArrivalTime = mean(Owls$ArrivalTime) )

calls.pred2 <- predict(m.nb2, New.data2, re.form=NA, SE.fit = TRUE)
# Error in predict.glmmTMB(m.nb2, New.data2, re.form = NA, SE.fit = TRUE) : 
#   re.form not yet implemented


 
##Compare to lme4 fit:
library(lme4)

nb4<-glmer.nb(NCalls~FoodTreatment+ArrivalTime+ 
                +(1|Nest), 
              data=Owls)

#Convergence warning
nb4.lsm<-lsmeans(nb4, ~FoodTreatment)

plot(nb4.lsm) #well-separated, glmmADMB SE's seem to be anomalous?



######### Follow-up: ##########################

##Check Russell Lenth's function to extract data from glmmADMB
#http://permalink.gmane.org/gmane.comp.lang.r.lme4.devel/12053
#Similar to Ben Bolker's function to extract from glmmTMB
#https://github.com/glmmTMB/glmmTMB/issues/205

#I cannot figure out why they give comparable results for pairwise comparisons,
  #but different results for standard errors
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

################################################################
#######Options to test significance of individual terms ###########
################################################################

#Car package
#https://r-forge.r-project.org/sendmessage.php?touser=1105

library(car)
Anova(m.nb2)
# Error in I.p[c(subs.relatives, subs.term), , drop = FALSE] : 
#   subscript out of bounds
# In addition: Warning message:
#   In is.na(coef(mod)) :
#   is.na() applied to non-(list or vector) of type 'NULL'

class(m.nb2)
class(m.nb)

nb2.reclass<-m.nb2
class(nb2.reclass)<-"glmmadmb"

class(nb2.reclass)

Anova(nb2.reclass)
#Error in t.default(as.vector(Y)) : argument is not a matrix

class(nb2.reclass)<-'merMod'
Anova(nb2.reclass)
# Error in vcov.merMod(mod) : 
#   trying to get slot "optinfo" from an object (class "merMod") that is not an S4 object 
?drop1
# Compute all the single terms in the scope argument that can be added to or dropped from the model, 
#fit those models and compute a table of the changes in fit.
drop1(m.nb2, test = "Chisq")
# Df    AIC    LRT  Pr(>Chi)    
# <none>           3478.4                     
# FoodTreatment  1 3516.2 39.795 2.821e-10 ***
#   ArrivalTime    1 3496.4 19.984 7.809e-06 ***

Anova(m.nb)
# Df  Chisq Pr(>Chisq)    
# FoodTreatment  1 41.935  9.437e-11 ***
#   ArrivalTime    1 20.878  4.895e-06 ***
#slightly different from drop1()


drop1(m.nb, test = "Chisq") #takes forever, refits model
# Df    AIC   LRT  Pr(>Chi)    
# <none>           3478.4                    
# FoodTreatment  1 3516.2 39.80 2.813e-10 ***
#   ArrivalTime    1 3496.4 19.98 7.826e-06 ***

#Same as drop1() for TMB model...
 #however Anova() and drop1() are not quite equivalent
 #Have to use drop1() if we use the TMB model

#Bootstrap with car?
?Boot
lm.simp<- lm(NCalls~FoodTreatment+ArrivalTime, 
             data=Owls)
Boot(lm.simp, R=500)
# original       bias    std. error
# t1* 29.3231729 -0.327651027   3.2481566
# t2* -3.2302850  0.009542173   0.5513673
# t3* -0.8522242  0.011749704   0.1274327
summary(lm.simp) #pretty close
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            29.3232     3.3589   8.730  < 2e-16 ***
#   FoodTreatmentSatiated  -3.2303     0.5158  -6.263 7.23e-10 ***
#   ArrivalTime            -0.8522     0.1345  -6.335 4.69e-10 ***