## ----setup---------------------------------------------------------------
require(GGally)

## ----data----------------------------------------------------------------
library(MASS)
data(UScrime)

## ----transform-----------------------------------------------------------
UScrime[,-2] = log(UScrime[,-2])

## ----bas-----------------------------------------------------------------
library(BAS)
crime.ZS =  bas.lm(y ~ ., 
                   data=UScrime,
                   prior="ZS-null",
                   modelprior=uniform(), initprobs="eplogp") 

## ---- fig.show='hold'----------------------------------------------------
plot(crime.ZS, ask=F)


## ----pip, fig.width=5, fig.height=5--------------------------------------
plot(crime.ZS, which = 4, ask=FALSE, caption="", sub.caption="")

## ----print---------------------------------------------------------------
crime.ZS

## ----summary-------------------------------------------------------------
summary(crime.ZS)

## ----image, fig.width=5, fig.height=5------------------------------------
image(crime.ZS, rotate=F)

## ----plot----------------------------------------------------------------
coef.ZS = coef(crime.ZS)
plot(coef.ZS, subset=c(5:6),  ask=F)

## ----coefall-------------------------------------------------------------

plot(coef.ZS, ask=FALSE)

## ----confint-coef--------------------------------------------------------

confint(coef.ZS)

## ----choice of estimator-------------------------------------------------
muhat.BMA = fitted(crime.ZS, estimator="BMA")
BMA  = predict(crime.ZS, estimator="BMA")

# predict has additional slots for fitted values under BMA, predictions under each model
names(BMA)

## ---- fig.width=5, fig.height=5------------------------------------------
par(mar=c(9, 9, 3, 3))
plot(muhat.BMA, BMA$fit, 
     pch=16, 
     xlab=expression(hat(mu[i])), ylab=expression(hat(Y[i])))
abline(0,1)

## ----HPM-----------------------------------------------------------------
HPM = predict(crime.ZS, estimator="HPM")

# show the indices of variables in the best model where 0 is the intercept
HPM$bestmodel

## ------------------------------------------------------------------------
(crime.ZS$namesx[HPM$bestmodel +1])[-1]

## ----MPM-----------------------------------------------------------------
MPM = predict(crime.ZS, estimator="MPM")
(crime.ZS$namesx[attr(MPM$fit, 'model') +1])[-1]

## ----BPM-----------------------------------------------------------------
BPM = predict(crime.ZS, estimator="BPM")
(crime.ZS$namesx[attr(BPM$fit, 'model') +1])[-1]

## ------------------------------------------------------------------------
library(GGally)
ggpairs(data.frame(HPM = as.vector(HPM$fit),  #this used predict so we need to extract fitted values
                   MPM = as.vector(MPM$fit),  # this used fitted
                   BPM = as.vector(BPM$fit),  # this used fitted
                   BMA = as.vector(BMA$fit))) # this used predict

## ----se------------------------------------------------------------------
BPM = predict(crime.ZS, estimator="BPM", se.fit=TRUE)
crime.conf.fit = confint(BPM, parm="mean")
crime.conf.pred = confint(BPM, parm="pred")
cbind(BPM$fit, crime.conf.fit, crime.conf.pred)


## ----MCMC----------------------------------------------------------------
crime.ZS =  bas.lm(y ~ ., 
                   data=UScrime,
                   prior="ZS-null",
                   modelprior=uniform(),
                   method = "MCMC") 

## ----diagnostics---------------------------------------------------------
diagnostics(crime.ZS, type="pip",  pch=16)
diagnostics(crime.ZS, type="model",  pch=16)

## ----biggerMCMC----------------------------------------------------------
crime.ZS =  bas.lm(y ~ ., 
                   data=UScrime,
                   prior="ZS-null",
                   modelprior=uniform(),
                   method = "MCMC", MCMC.iterations = 10^6)  

diagnostics(crime.ZS, type="model", pch=16)

