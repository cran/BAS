% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EB_global.R
\name{EB.global}
\alias{EB.global}
\alias{EB.global.bas}
\title{Find the global Empirical Bayes estimates for BMA}
\usage{
EB.global(object, tol = 0.1, g.0 = NULL, max.iterations = 100)
}
\arguments{
\item{object}{A 'bas' object created by \code{\link{bas}}}

\item{tol}{tolerance for estimating g}

\item{g.0}{initial value for g}

\item{max.iterations}{Maximum number of iterations for the EM algorithm}
}
\value{
An object of class 'bas' using Zellner's g prior with an estimate of
g based on all models
}
\description{
Finds the global Empirical Bayes estimates of g in Zellner's g-prior and
model probabilities
}
\details{
Uses the EM algorithm in Liang et al to estimate the type II MLE of g in
Zellner's g prior
}
\examples{

library(MASS)
data(UScrime)
UScrime[,-2] = log(UScrime[,-2])
# EB local uses a different g within each model
crime.EBL =  bas.lm(y ~ ., data=UScrime, n.models=2^15,
                    prior="EB-local", initprobs= "eplogp")
# use a common (global) estimate of g
crime.EBG = EB.global(crime.EBL)


}
\references{
Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger, J.O.
(2008) Mixtures of g-priors for Bayesian Variable Selection. Journal of the
American Statistical Association. 103:410-423.  \cr
\doi{10.1198/016214507000001337}
}
\seealso{
\code{\link{bas}}, \code{\link{update}}
}
\author{
Merlise Clyde \email{clyde@stat.duke.edu}
}
\concept{coef priors}
\keyword{regression}
