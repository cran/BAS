% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confint.R
\name{confint.pred.bas}
\alias{confint.pred.bas}
\title{Compute Credible (Bayesian Confidence) Intervals for a BAS predict object}
\usage{
\method{confint}{pred.bas}(object, parm, level = 0.95, nsim = 10000, ...)
}
\arguments{
\item{object}{an object created by \code{\link{predict.bas}}}

\item{parm}{character variable, "mean" or "pred".  If missing parm='pred'.}

\item{level}{the nominal level of the (point-wise) credible interval}

\item{nsim}{number of Monte Carlo simulations for sampling methods with BMA}

\item{...}{optional arguments to pass on to next function call; none at this
time.}
}
\value{
a matrix with lower and upper level * 100 percent credible intervals
for either the mean of the regression function or predicted values.  %%
}
\description{
Compute credible intervals for in-sample or out of sample prediction or for
the regression function
}
\details{
This constructs approximate 95 percent Highest Posterior Density intervals
for 'pred.bas' objects.  If the estimator is based on model selection, the
intervals use a Student t distribution using the estimate of g.  If the
estimator is based on BMA, then nsim draws from the mixture of Student t
distributions are obtained with the HPD interval obtained from the Monte
Carlo draws. %% ~~ If necessary, more details than the description above ~~
}
\examples{

data("Hald")
hald.gprior =  bas.lm(Y~ ., data=Hald, alpha=13, prior="g-prior")
hald.pred = predict(hald.gprior, estimator="BPM", predict=FALSE, se.fit=TRUE)
confint(hald.pred, parm="mean")
confint(hald.pred)  #default
hald.pred = predict(hald.gprior, estimator="BMA", predict=FALSE, se.fit=TRUE)
confint(hald.pred)


}
\seealso{
\code{\link{predict.bas}}

Other bas methods: 
\code{\link{BAS}},
\code{\link{bas.lm}()},
\code{\link{coef.bas}()},
\code{\link{confint.coef.bas}()},
\code{\link{diagnostics}()},
\code{\link{fitted.bas}()},
\code{\link{force.heredity.bas}()},
\code{\link{image.bas}()},
\code{\link{plot.confint.bas}()},
\code{\link{predict.bas}()},
\code{\link{predict.basglm}()},
\code{\link{summary.bas}()},
\code{\link{update.bas}()},
\code{\link{variable.names.pred.bas}()}

Other CI methods: 
\code{\link{confint.coef.bas}()},
\code{\link{plot.confint.bas}()}
}
\author{
Merlise A Clyde
}
\concept{CI methods}
\concept{bas methods}
\keyword{regression}
