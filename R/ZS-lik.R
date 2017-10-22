g.lik = function(g, R2, n, p, log=TRUE) {
 out = .5*(log(1 + g)*(n - p - 1.0) - (n-1)*log(1 + (1 - R2)*g))
 if (any(is.na(out))) browser()
 if (! log) out=exp(out)
 return(out)
  }

ZS.lik.g = function(g, R2, n, p, rscale=1, log=TRUE) {
  out=.5*(log(.5*n*rscale) -3*log(g) - rscale*n/g) - log(gamma(.5)) +
   g.lik(g, R2, n, p, log=TRUE)
  if (!log) out=exp(out)
  return(out)
}

ZS.lik.tau = function(tau, R2, n, p, rscale=1, log=TRUE) {
  out= #.5*(log(n*rscale*.5) - log(tau) - rscale*n*tau) - log(gamma(.5)) +
     dgamma(tau, shape=.5, rate=n/2, log=TRUE) +
      g.lik(1/tau, R2, n, p, log=TRUE)
  if (!log) out=exp(out)
  return(out)
}
ZS.marglik = function(R2, n, p, rscale=1, lower=0, upper=Inf, invgamma=TRUE) {
  if (invgamma)   integrate(f=ZS.lik.tau, lower, upper, R2, n, p, rscale=1, log=FALSE)
  else   integrate(f=ZS.lik.g, lower=0, upper=Inf, R2, n, p, rscale=1, log=FALSE)
}
