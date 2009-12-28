bas.lm = function(formula, data, n.models=NULL,  prior="ZS-null", alpha=NULL,
                  modelprior=uniform(),
                  initprobs="Uniform", random=TRUE, method="BAS", update=NULL, 
                  bestmodel=NULL, bestmarg=NULL, prob.local=0.0)  {
  call = match.call()
  lm.obj = lm(formula, data, y=TRUE, x=TRUE)
  Y = lm.obj$y
  X = lm.obj$x
  p = dim(X)[2]

  
  if (!is.numeric(initprobs)) {
    initprobs = switch(initprobs,
      "Uniform"= c(1.0, rep(.5, p-1)),
      "eplogp" = eplogprob(lm.obj)
      )
  }
   if (length(initprobs) == (p-1))
     initprobs = c(1.0, initprobs)
   if (length(initprobs) != p)
    stop(simpleError(paste("length of initprobs is not", p)))

  pval = summary(lm.obj)$coefficients[,4]
  if (any(is.na(pval))) {
    print(paste("full model is rank deficient; automatically dropping",
                  sum(is.na(pval)),
                  " variables."))
    initprobs[is.na(pval)] = 0.0
  }

  if (initprobs[1] < 1.0) initprobs[1] = 1.0
# intercept is always included otherwise we get a segmentation
# fault (relax later)

  prob = as.numeric(initprobs)

  if (is.null(n.models)) n.models = 2^(p-1)
  deg = sum(initprobs >= 1) + sum(initprobs <= 0)
  if (deg > 1 & n.models == 2^(p - 1)) {
    n.models = 2^(p - deg)
    print(paste("There are", as.character(deg),
                "degerate sampling probabilities (0 or 1); decreasing the number of models to",                 as.character(n.models)))
  }

  if (n.models > 2^30) stop("Dimension of model space is too big to enumerate\n  Rerun with a smaller value for n.models")
  if (n.models > 2^20)
    print("Number of models is REALLY BIG -this may take a while")


  if (modelprior$family == "Bernoulli") {
   if (length(modelprior$hyper.parameters) == 1) 
      modelprior$hyper.parameters = c(1, rep(modelprior$hyper.parameters, p-1))
    if  (length(modelprior$hyper.parameters) == (p-1)) 
     modelprior$hyper.parameters = c(1, modelprior$hyper.parameters)
    if  (length(modelprior$hyper.parameters) != p)
      stop(" Number of probabilities in Bernoulli family is not equal to the number of variables or 1")
  }
  
  int = TRUE  # assume that an intercept is always included 
  method.num = switch(prior,
    "g-prior"=0,
    "hyper-g"=1,
    "EB-local"=2,
    "BIC"=3,
    "ZS-null"=4,
    "ZS-full"=5,
    "hyper-g-laplace"=6,
    "AIC"=7,
    "EB-global"=2,
    "hyper-g-n"=8
    )
  if (is.null(alpha) &&
      (method.num == 0 || method.num == 1 || method.num  == 6)) {
    stop(simpleError(paste("Must specify a value of alpha for", prior)))
  }

  if (is.null(alpha)) alpha=0.0 
  if (is.null(bestmodel)) {
    bestmodel = as.integer(initprobs)
    bestmarg = -Inf}
  if (is.null(bestmarg)) bestmarg = 0
  if (is.null(update)) update = n.models+1

  modelindex = as.list(1:n.models)
  Yvec = as.numeric(Y)
  R2 = as.numeric(rep(0.0, n.models))
  beta = as.list(1:n.models)
  se = as.list(1:n.models)
  mse = as.numeric(rep(0.0, n.models))
  modelspace = as.list(1:n.models)
  modelprobs = as.numeric(rep(0.0, n.models))
  priorprobs = as.numeric(rep(1.0, n.models))
  logmargy = as.numeric(rep(0.0, n.models))
  shrinkage = as.numeric(rep(0.0, n.models))
  modeldim = as.integer(rep(0, n.models))
  sampleprobs = as.double(rep(0.0, n.models))
  modeltree = list(NULL, NULL, NULL, NULL, FALSE)
  if (random) { 
  if (method == "BAS")
     ANS = .Call("sampleworep",
      Yvec, X,
      prob, R2,beta, se, mse, modelspace, modelprobs,
      priorprobs,logmargy, sampleprobs,
      modeldim, shrinkage, incint=as.integer(int), 
      alpha= as.numeric(alpha),
      method=as.integer(method.num), modelprior=modelprior,
      update=as.integer(update),
      Rbestmodel=as.integer(bestmodel),
      Rbestmarg=as.numeric(bestmarg),
      plocal=as.numeric(prob.local),
      PACKAGE="BAS")
      prob = ANS[[1]]
}
  else {
    ans = .Call("deterministic",
      Yvec, X,
      prob,
      R2,beta, se, mse, modelspace,  modelprobs,
      priorprobs,logmargy, sampleprobs,
      modeldim, shrinkage, incint=as.integer(int),
      alpha= as.numeric(alpha),
      method=as.integer(method.num),modelprior=modelprior,
      PACKAGE="BAS")
  }

  namesx = dimnames(X)[[2]]
  namesx[1] = "Intercept"
  n.models = as.integer(n.models)
  result =  list(namesx=namesx, which=modelspace, postprobs=modelprobs, priorprobs=priorprobs,
    logmarg=logmargy, R2=R2, mse=mse, ols=beta, ols.se=se,
    shrinkage=shrinkage, probne0=prob, size=modeldim, n=length(Yvec),
    prior=prior, modelprior=modelprior,alpha=alpha, n.models=n.models, n.vars=p,
    sampleprobs=sampleprobs, Y=Yvec, X=X, call=call)

  class(result) = "bma"

  if (prior == "EB-global") result = EB.global.bma(result)
    
  return(result) 
}
