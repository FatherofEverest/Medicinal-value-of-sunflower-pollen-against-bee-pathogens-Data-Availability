######Try to fit a model that can deal with complete separation #####
## This could be used to deal with pupal production

#install.packages("blme") #Bayesian linear mixed-effects...
# library(blme)
# cmod_blme_L2 <- bglmer(predation~ttt+(1|block),data=newdat,
#                        family=binomial,
#                        fixef.prior = normal(cov = diag(9,4)))
# There are other packages in R (brglm, logistf) 
# that can handle completely separated data, 
# but they only apply to logistic regression, 
# and they cannot simultaneously incorporate random effects 
# (Pasch et al., American Naturalist 2013 
#   used brglm to handle completely separated data, arguing that random effects were not very important in their system.)

library(blme)
Data <- data.frame (Response = c(rep(1, 6), rep(0, 14), rep(0, 20)),
                                 Treatment = c(rep("Gooddiet", 20), rep("Baddiet", 20)),
                                 Colony = c("C1", "C2"))
Data
?bglmer #Maximum a posteriori estimation for 
       #linear and generalized linear mixed-effects models in a Bayesian setting. 
          ##Built off of lmer.
bmod1 <- bglmer(Response~Treatment+(1|Colony),data=Data, family = binomial,
                fixef.prior=normal)
                      
summary(bmod1)
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)         -3.820      1.494  -2.556   0.0106 *
#   TreatmentGooddiet    2.812      1.305   2.155   0.0312 *

#Try significance tests:

library(car)
Anova(bmod1)
# Chisq Df Pr(>Chisq)  
# Treatment 4.6429  1    0.03118 *

drop1(bmod1) #nope

library(lsmeans)
lsm.bmod1<-lsmeans(bmod1, ~Treatment) #still works!
plot(lsm.bmod1)

#back-transform
library(boot)
lsm.df<- as.data.frame(summary(lsm.bmod1))
lsm.df$btmean <-inv.logit(lsm.df$lsmean)
lsm.df$bthi<- inv.logit(lsm.df$lsmean + lsm.df$SE)
lsm.df$btlo<- inv.logit(lsm.df$lsmean - lsm.df$SE)

library(ggplot2)
ggplot(lsm.df, aes(x=Treatment, y=btmean, shape = Treatment))+
  geom_pointrange(aes(y=btmean, ymax=bthi, ymin=btlo))
