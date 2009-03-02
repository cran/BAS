\name{bas.lm}
\alias{bas}
\alias{bas.lm}
\title{Bayesian Adaptive Sampling Without Replacement for Variable
  Selection in Linear Models}
\description{Sample without replacement from a posterior distribution on models}
\usage{bas.lm(formula, data, n.models, alpha=NULL, 
 prior="ZS", initprobs="Uniform", random=TRUE, update=NULL,
 bestmodel = NULL, bestmarg = NULL, prob.local = 0, user.prob=NULL) }
\arguments{
  \item{formula}{linear model formula for the full model with all
    predictors, Y ~ X.  All code assumes that an intercept will be
    included in each model.}
  \item{data}{data frame}
  \item{n.models}{number of models to sample}
  \item{initprobs}{vector of initial marginal inclusion
    probabilities used for sampling without replacement or method,
    if "Uniform" each predictor variable is equally likely to be
    sampled (equivalent to random sampling without replacement). If
    "eplogp", use the \code{\link{eplogprob}} function to aproximate the 
    Bayes factor to find initial marginal inclusion probabilitites and
    sample without replacement the model probabilities using these
    inclusion probabilaties. If "user", the user may specify a vector of
  probabilities using the optional argument user.prob.  This should
  include the intercept with probability one.} 
  \item{alpha}{optional hyperparameter in g-prior or hyper g-prior.  For Zellner's
    g-prior, alpha = g, for the Liang et al hyper-g method, recommended choice is alpha
    = 3 or 4.}
  \item{prior}{prior distribution for regression coefficients.  Choices include
"AIC", "BIC", "g-prior", "ZS-null", "ZS-full", "hyper-g", "hyper-g-laplace", "EB-local", and "EB-global" }
  \item{random}{A logical variable indicating whether to use the
    stochastic (random=TRUE) or deterministic (random=FALSE) algorithm
    for sampling models without replacement}
  \item{update}{number of iterations between potential updates of the
    sampling probabilities. If NULL do not update, otherwise the
    algorithm will update using the marginal inclusion probabilities}
  \item{bestmodel}{optional binary vector representing a model to
    initialize the sampling. If  NULL sampling starts with the Full
    model}
  \item{bestmarg}{optional value for the log marginal associated with
    the bestmodel}
  \item{prob.local}{An experimental option to allow sampling of models
    "near" the median probability model.  Not recommended for use at
    this time}
  \item{user.prob}{User specified probabilities for sampling}
}
\details{
BAS provides two search algorithms to find high probability models for
use in  Bayesian Model Averaging or Bayesian model selection. 
}
\value{
  \code{bas} returns an object of class \code{BMA}

An object of class \code{BMA} is a list containing at least the following components:

  \item{postprob}{the posterior probabilities of the models selected}
  \item{namesx}{the names of the variables}
  \item{R2}{R2 values for the models}
  \item{logmarg}{values of the log of the marginal likelihood for the
    models}
 \item{n.vars}{total number of independent variables in the full model,
   including the intercept}
  \item{size}{the number of independent variables in each of the models,
    includes the intercept}
  \item{which}{a list of lists with one list per model with  variables
    that are included in the model} 
  \item{probne0}{the posterior probability that each variable is non-zero}
  \item{ols}{list of lists with one list per model giving the OLS
estimate of each (nonzero) coefficient for each model}
  \item{ols.se}{list of lists with one list per model giving the OLS standard error of each coefficient for each model}
  \item{prior}{the name of the prior  that created the BMA object}
  \item{alpha}{value of hyperparameter in prior used to create the BMA
object.}
  \item{Y}{response}
  \item{X}{matrix of predictors}
  
 The function  \code{\link{summary.bma}}, is used to print a summary of
 the results. The function \code{\link{plot.bma}} is used to plot
 posterior distributions for the coefficients and
 \code{\link{image.bma}} provides an image of the distribution over models.
 Posterior summaries of coefficients can be extracted using
 \code{\link{coefficients.bma}}.  Fitted values and predictions can be
 obtained using the functions  \code{\link{fitted.bma}} and \code{\link{predict.bma}}.
 BMA objects may be updated to use a different prior (without rerunning
 the sampler) using the function  \code{\link{update.bma}}.
}

\references{
  Clyde, M. Ghosh, J. and Littman, M. (2009) Bayesian Adaptive Sampling
  for Variable Selection. Department of Statistical Science Discussion
  Paper. Duke University.
  
  Clyde, M. and George, E. I. (2004) Model Uncertainty. Statist. Sci.,
  19, 81-94. \cr \url{http://www.isds.duke.edu/~clyde/papers/statsci.pdf}
  
  Clyde, M. (1999)
  Bayesian Model Averaging and Model Search Strategies (with
  discussion). In Bayesian Statistics 6. J.M. Bernardo, A.P. Dawid,
  J.O. Berger, and A.F.M. Smith eds. Oxford University Press, pages
  157-185.

  Hoeting, J. A., Madigan, D., Raftery, A. E. and Volinsky, C. T. (1999)
  Bayesian model averaging: a tutorial (with discussion). Statist. Sci.,
  14, 382-401. \cr \url{http://www.stat.washington.edu/www/research/online/hoeting1999.pdf}
  
  Liang, F., Paulo, R., Molina, G., Clyde, M. and  Berger,
  J.O. (2005) Mixtures of  g-priors for Bayesian Variable
  Selection.  \cr \url{http://www.stat.duke.edu/05-12.pdf}
  
  Zellner, A. (1986) On assessing prior distributions and Bayesian
  regression analysis with g-prior distributions. In Bayesian Inference
  and Decision Techniques: Essays in Honor of Bruno de Finetti,
  pp. 233-243. North-Holland/Elsevier.
    
  Zellner, A. and Siow, A. (1980) Posterior odds ratios for selected
  regression hypotheses. In Bayesian Statistics: Proceedings of the First
  International Meeting held in Valencia (Spain), pp. 585-603.
}


\author{Merlise Clyde (\email{clyde@stat.duke.edu}) and Michael Littman}
  
\seealso{ \code{\link{summary.bma}},
          \code{\link{coefficients.bma}},
          \code{\link{print.bma}},
	  \code{\link{predict.bma}},
	  \code{\link{fitted.bma}}
          \code{\link{plot.bma}},
	  \code{\link{image.bma}},
	  \code{\link{eplogprob}},
	  \code{\link{update.bma}}
	  
}
	
\examples{
demo(BAS.hald)
\dontrun{demo(BAS.USCrime) }
}

\note{Uniform prior probabilities on  models are the only option
  currently.  A future update should allow alternative priors on models
  to be incorporated into the sampling and posterior inference.  For
  now, users may manually reweight output using the log marginal
  likelihoods to update posterior model probabilities and probne0.
}
\keyword{regression}
\keyword{models}

