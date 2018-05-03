
# BAS 1.5.0  May 2, 2018

## Features

* add na.action for handling NA's for predict methods 
[issue #10]( https://github.com/merliseclyde/BAS/issues/10)

* added `include.always` as new argument to `bas.lm`.  This allows a formula to specify which terms should always be included in all models.  By default the intercept is always included.

* added a section to the vignetted to illustrate weighted regression and the 
`force.heredity.bas` function to group levels of a factor so that they enter 
or leave the model together.



## Bugs

* fixed problem if there is only one model for `image` function;  
github [issue #11](https://github.com/merliseclyde/BAS/issues/11)

* fixed error in `bas.lm` with non-equal weights where R2 was incorrect.
 [issue #17](https://github.com/merliseclyde/BAS/issues/17)
## Deprecated 
 
* deprecate the `predict` argument in `predict.bas`, `predict.basglm` and internal functions as it is not utilized


#  BAS 1.4.9 March 24, 2018

## Bugs

* fixed bug in `confint.coef.bas` when parm is a character string
* added parentheses in betafamily.c line 382 as indicated in CRAN check for R devel
 
## Features

* added option to determine k for `Bayes.outlier` if prior 
probability of no outliers is provided

#  BAS 1.4.8 March 10, 2018

## Bugs 

* fixed issue with scoping in eval of data in `predict.bas` if dataname is defined in local env.

* fixed issue 10 in github (predict for estimator='BPM'
failed  if there were NA's in the X data.  Delete NA's
before finding the closest model.

* fixed bug in 'JZS' prior - merged pull request #12 from vandenman/master

* fixed bug in bas.glm when default betaprior (CCH) is used and inputs were INTEGER instead of REAL

* removed warning with use of 'ZS-null' for backwards compatibility

## Features added

* updated print.bas to reflect changes in print.lm

* Added Bayes.outlier function to calculate posterior probabilities of outliers using the method from 
Chaloner & Brant for linear models.

#  BAS 1.4.7 October 22, 2017

## Updates

* Added new method for `bas.lm` to obtain marginal likelihoods with the Zellner-Siow Priors for "prior= 'JZS' using QUADMATH routines for numerical integration.  The optional hyperparameter alpha may now be used to adjust the scaling of the ZS prior where g ~ G(1/2, alpha*n/2) as in the `BayesFactor` package of Morey, with a default of alpha=1 corresponding to the ZS prior used in Liang et al (2008).  This also uses more stable evaluations of log(1 + x) to prevent underflow/overflow.

* Priors `ZS-full` for bas.lm is planned to be deprecated.  

* replaced math functions to use portable C code from Rmath and consolidated header files

#  BAS 1.4.6 May 24, 2017

## Updates

*  Added force.heredity.interaction function to allow higher order interactions to be included only if their "parents" or lower order interactions or main effects were included.   Currently tested with two way interactions.  This is implemented post-sampling; future  updates will add this at the sampling stage which will reduce memory usage and sampling times by reducing the number of models under consideration.

## Bugs

* Fixed unprotected ANS in C code in glm_sampleworep.c and sampleworep.c after call to PutRNGstate and possible stack imbalance in glm_mcmc.

* Fixed problem with predict for estimator=BPM when newdata was one row 


#  BAS 1.4.5 March 28, 2017

## Bugs

* Fixed non-conformable error with `predict` when new data was from a dataframe with one row.

* Fixed problem with missing weights for prediction using the median probability model with no new data.

#  BAS 1.4.4 March 14, 2017

## Updates 

* Extract coefficent summaries, credible intervals and plots for the `HPM` and ` MPM` in addition to the default `BMA` by adding a new `estimator` argument to the `coef` function. The new `n.models` argument to `coef` provides summaries based on the top `n.models` highest probability models to reduce computation time. 'n.models = 1' is equivalent to the highest probability model.

* use of newdata that is a vector is now depricated for predict.bas; newdata must be a dataframe or missing, in which case fitted values based on the dataframe used in fitting is used

* factor levels are handled as in `lm` or `glm` for prediction when there may be only level of a factor in the newdata

## Bugs

* fixed issue for prediction when newdata has just one row

* fixed missing id in plot.bas for which=3

#  BAS 1.4.3  February 18, 2017

## Updates

* Register symbols for foreign function calls
* bin2int is now deprecated
* fixed default MCMC.iteration in `bas.lm` to agree with documentation
* updated vignette to include more examples, outlier detection, and finding the best predictive probability model
* set a flag for MCMC sampling `renormalize` that selects whether the Monte Carlo frequencies are used to estimate posterior model and marginal inclusion probabilities (default `renormalize = FALSE`) or that marginal likelihoods time prior probabilities that are renormalized to sum to 1 are used.  (the latter is the only option for the other methods); new slots for probne0.MCMC, probne0.RN, postprobs.RN and postprobs.MCMC.

## Bug fixes

 *  fixed problem with prior.bic, robust, and hyper.g.n where default had missing n that was not set in hyperparameters
 * fixed error in predict and plot for GLMs when family is provided as a function
 

# BAS 1.4.2   October 12, 2016

## Updates

* added df to the object returned by bas.glm to simplify `coefficients` function.

## Bug Fixes
* corrected expected value of shrinkage for intrinsic, hyper-g/n and TCCH priors for glms

# BAS 1.4.1   September 17, 2016

## Bug Fixes

* the modification in 1.4.0 to automatically handle NA's led to
  errors if the response was transformed as part of the forumula; 
  this is fixed 

## Features
* added subset argument to `bas.lm` and `bas.glm`

# BAS 1.4.0   August 25, 2016

## New features

* added `na.action` for `bas.lm` and `bas.glm` to omit missing data.
* new function to plot credible intervals created by `confint.pred.bas` or `confint.coef.bas`.   See the help files for an example or the vignette.
* added `se.fit` option in `predict.basglm`.
* Added `testBF` as a `betaprior` option for `bas.glm` to implement Bayes Fatcors based on the likelihood ratio statistic's distribution for GLMs.
* DOI for this version is http://dx.doi.org/10.5281/zenodo.60948


# BAS 1.3.0   July 15, 2016

## New Features

A vignette has been added at long last!  This illustrates several of the new features in `BAS` such as
  
* new functions for computing credible intervals for fitted and predicted values  `confint.pred.bas()`
* new function for adding credible intervals for coefficients  `confint.coef.bas()` 
* added posterior standard deviations for fitted values and predicted values in  `predict.bas()`
  


## Deprication 
* deprecated use of `type` to specify estimator in fitted.bas	and replaced with `estimator` so that `predict()` and `fitted()` are compatible with other S3 methods. 
* updated funtions to be of class `bas` to avoid NAMESPACE conficts with  other libraries
	

# BAS 1.2.2 June 29, 2016

## New Features
* added option to find "Best Predictive Model" or "BPM" for `fitted.bas` or `predict.bas`
* added local Empirical Bayes prior and fixed g-prior for `bas.glm`
* added `diagnostic()` function for checking convergence of `bas` objects created with `method = "MCMC"`"
* added truncated power prior as in Yang, Wainwright & Jordan (2016)

##	Minor Changes

* bug fix in `plot.bas` that appears with Sweave
* bug fix in `coef.bma` when there is just one predictor
	
	

# BAS 1.2.1  April 16, 2016
* bug fix for method="MCMC" with truncated prior distributions
	where MH ratio was incorrect allowing models with 0 probability to
	be sampled.
* fixed error in Zellner-Siow prior (ZS-null) when n=p+1 or
	saturated model  where log marginal likelihood should be 0 

# BAS 1.2.0  April 11, 2016
* removed unsafe code where Rbestmarg (input) was being
	overwritten in .Call which would end up in corruption of the
	constant pool of the byte-code  (Thanks to Tomas Kalibera for
	catching this!)
* fixed issue with dimensions for use with Simple Linear Regression

# BAS 1.1.0   March 31, 2016

## New Features
* added truncated Beta-Binomial prior and truncated Poisson (works
	only with MCMC currently)
* improved code for finding fitted values under the Median 
* deprecated method = "AMCMC" and issue warning message

## Minor Changes
* Changed S3 method for plot and image to use class `bas` rather than
	 `bma` to avoid name conflicts with other packages

# BAS 1.09
	- added weights for linear models
	- switched LINPACK calls in bayesreg to LAPACK finally should be
	faster
	- fixed bug in intercept calculation for glms
	- fixed inclusion probabilities to be a vector in the global EB
	methods for linear models
# BAS 1.08
	- added intrinsic prior for GLMs
	- fixed problems for linear models for p > n and R2 not correct
# BAS 1.07
	- added phi1 function from Gordy (1998)  confluent hypergeometric
	function of two variables  also known as one of the Horn
	hypergeometric functions or Humbert's phi1
	- added Jeffrey's prior on g
	- added the general tCCH prior and special cases of the hyper-g/n.
	- TODO check shrinkage functions for all	
# BAS 1.06
	- new improved Laplace approximation for hypergeometric1F1
	- added class basglm for predict
	- predict function now handles glm output
	- added dataframe option for newdata in predict.bas and predict.basglm
	- renamed coefficients in output to be 'mle' in bas.lm to be consistent across
	lm and glm versions so that predict methods can handle both
	cases.  (This may lead to errors in other external code that
	expects object$ols or object$coefficients)
	- fixed bug with initprobs that did not include an intercept for bas.lm
	
# BAS 1.05
	- added thinning option for MCMC method for bas.lm
	- returned posterior expected shrinkage for bas.glm
	- added option for initprobs = "marg-eplogp" for using marginal
	SLR models to create starting probabilities or order variables
	especially for p > n case
	- added standalone function for hypergeometric1F1 using Cephes
	library and a Laplace aproximation
	-Added class "BAS" so that predict and fitted functions (S3
	methods) are not masked by functions in the BVS package: to do
	modify the rest of the S3 methods.
	
# BAS 1.04

	- added bas.glm for model averaging/section using mixture of g-priors for
	GLMs.  Currently limited to Logistic Regression
	- added Poisson family for glm.fit

# BAS 1.0	
	- cleaned up  MCMC method code
	
# BAS 0.93

	- removed internal print statements in bayesglm.c
	- Bug fixes in AMCMC algorithm

# BAS 0.92

	- fixed glm-fit.R  so that hyperparameter for BIC is numeric

# BAS 0.91

	- added new AMCMC algorithm

# BAS 0.91

	- bug fix in bayes.glm

# BAS 0.90

	- added C routines for fitting glms

# BAS 0.85

	- fixed problem with duplicate models if n.models was > 2^(p-1) by
   restricting n.models

	- save original X as part of object so that fitted.bma gives the
   correct fitted values (broken in version 0.80)
 
# BAS 0.80

	- Added `hypergeometric2F1` function that is callable by R
	- centered X's in bas.lm so that the intercept has the correct
  shrinkage
	- changed `predict.bma` to center newdata using the mean(X)
	- Added new Adaptive MCMC option (method = "AMCMC")  (this is not
  stable at this point)

# BAS 0.7

	-Allowed pruning of model tree to eliminate rejected models
 
# BAS 0.6

	- Added MCMC option to create starting values for BAS (`method = "MCMC+BAS"`)

# BAS 0.5

	-Cleaned up all .Call routines so that all objects are duplicated or
 allocated within code

# BAS 0.45

	- fixed ch2inv that prevented building on Windows in bayes glm_fit

# BAS 0.4

	- fixed fortran calls to use F77_NAME macro 
	- changed  allocation of objects for .Call to prevent some objects from being overwritten.  
  
# BAS 0.3

	- fixed EB.global function to include prior probabilities on models
	- fixed update function 

# BAS 0.2

	- fixed predict.bma to allow newdata to be a matrix or vector with the
  column of ones for the intercept optionally included.
	- fixed help file for predict 
	- added modelprior argument to bas.lm so that users may now use the
	beta-binomial prior distribution on model size in addition to the
	default uniform distribution
	- added functions uniform(), beta-binomial() and Bernoulli() to create
	model prior objects
	- added a vector of user specified initial probabilities as an option for
	argument initprobs in bas.lm and removed the separate argument user.prob

