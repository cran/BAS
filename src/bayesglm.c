#include "sampling.h"

/* Version of glm.fit that can be called directly from R or C*/
SEXP glm_fit(SEXP RX, SEXP RY, SEXP Reta, SEXP Rcoef, 
	     SEXP Rse, SEXP Rdeviance, 
             SEXP Rweights, SEXP Rresiduals,
	     SEXP Reffects,
             SEXP Rpivot,  SEXP Rqrauxmat, SEXP Rworkmat)
	     

{
  SEXP RYwork, RXwork, Rstart;
  double  *Xwork, *Ywork, *coef,*start, *work, *qraux, *residuals,*effects,
           one, zero, tol;
  int  job, l, i, j, info, inc, n, p, rank, *pivot, *xdims;


  zero = 0.0;
  one = 1.0;
  inc = 1;
  job = 01;
  info = 1;
  xdims = INTEGER(getAttrib(RX,R_DimSymbol));
  n = xdims[0];
  p = xdims[1];
  rank = 1;
  tol = DBL_EPSILON;
 

  PROTECT(RYwork = coerceVector(RY, REALSXP));
  Ywork = REAL(RYwork);
  PROTECT(RXwork =  coerceVector(RX, REALSXP));
  Xwork = REAL(RXwork);
  PROTECT(Rstart = coerceVector(Rcoef, REALSXP));
  start = REAL(Rstart);

  pivot = INTEGER(Rpivot);
  coef = REAL(Rcoef);
  residuals = REAL(Rresiduals);
  qraux = REAL(Rqrauxmat);
  effects = REAL(Reffects);
  work = REAL(Rworkmat);
  
  /*  eta <- rep.int(0, nobs) + offset
  mu <- linkinv(eta);
  dev <- sum(dev.resids(y, mu, weights))
  w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
  residuals <- (y - mu)/mu.eta(eta)
  z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
  w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
  ngoodobs <- as.integer(nobs - sum(!good))
  */ 
  Rprintf("calling dqrls\n");
  dqrls_(&Xwork[0], &n, &p, &Ywork[0], &inc, &tol,  &start[0],
	 &residuals[0], &effects[0], &rank, &pivot[0], &qraux[0], &work[0]);
  for (j=0; j<p; j++) { coef[pivot[j]] = start[j];}
  /* eta = drop(x %*% start)+ offset;
  mu  = linkinv(eta);
  dev = sum(dev.resids(y, mu, weights));
  */
  UNPROTECT(2);
  
  return(Rcoef);
  //  dpofa_(&XtX[0],&p,&p, &info);
  //  dposl_(&XtX[0],&p,&p,&coefficients[0]);
  //  dpodi_(&XtX[0],&p,&p, &det, &job);

  //  ete = ddot_(&p, &coefficients[0], &inc, &XtY[0], &inc);
  //  *mse = (*mse - ete)/((double) (n - p)); 

  /*  for (j=0, l=0; j < p; j++)  {
    for (i=0; i <  p; i++) {
       if (i == j)  {
	se[j] = sqrt(XtX[l] * *mse);
       }
      l += 1;    }}
  */
}



/* Version of gexpectations that can be called directly from R */
void gexpectations_glm_vect(int *nmodels, int *p, int *pmodel, int *nobs, double *R2, double *alpha, int *method,
                        double *RSquareFull, double *SSY, double *logmarg, double *shrinkage) {

  int i; 
     for (i=0; i < *nmodels; i++) {
	 gexpectations(*p, pmodel[i],  *nobs,  R2[i], *alpha, *method, 
		        *RSquareFull,  *SSY,  &logmarg[i],  &shrinkage[i]);

	     }
}


void gexpectations_glm(int p, int pmodel, int nobs, double R2, double alpha, int method,  double RSquareFull, double SSY, double *logmarg, double *shrinkage) {  
    
    *shrinkage = 1.0;

    switch (method) { 
     case 0: 
      *logmarg = logBF_gprior(R2, nobs, pmodel,alpha); 
      if (pmodel > 1) *shrinkage = alpha/(1.0 + alpha);
      break;
     case 1:  
      *logmarg = logBF_hyperGprior(R2, nobs, pmodel,alpha); 
      *shrinkage = shrinkage_hyperg(R2, nobs, pmodel, alpha);
      break;
     case 2: 
      *logmarg = logBF_EB(R2, nobs, pmodel,alpha); 
      *shrinkage = shrinkage_EB_local(R2, nobs, pmodel,alpha);
      break;
     case 3: 
     *logmarg = BIC(R2, nobs, pmodel,SSY); 
      break;
     case 4: 
       *logmarg = LogBF_ZS_null(R2, nobs, pmodel);
       *shrinkage = E_ZS_approx_null(R2,nobs,pmodel-1);
      break; 
     case 5: 
       *logmarg = LogBF_ZS_full(RSquareFull, R2, nobs, p, pmodel);
      break; 
      case 6: 
        *logmarg =  logBF_hyperGprior_laplace(R2, nobs, pmodel,alpha); 
	*shrinkage = shrinkage_laplace(R2, nobs, pmodel, alpha); 
      break; 
      case 7: 
	*logmarg =  AIC(R2,  nobs, pmodel,SSY);
      break; 
      case 8: 
	  *logmarg =  LogBF_Hg_null(R2,  nobs, pmodel, alpha, 1);
 	  if (pmodel > 1) {
	      *shrinkage = LogBF_Hg_null(R2,  nobs, pmodel+2, alpha, 2);
	      *shrinkage = exp(*shrinkage - *logmarg); }
      break; 
  default: 
      Rprintf("Error: Method must be one of g-prior, hyper-g, laplace (hyper-g), AIC, BIC, ZS-null, or ZS-full\n");
      break;
    }
}






double logBF_gprior_glm( double Rsquare, int n,  int p, double g)
  {  double logmargy;
  /* log Marginal under reference prior for phi, intercept and
     g prior for regression coefficients 
     log marginal for null model is 0 */
  logmargy =  .5*(log(1.0 + g) * ((double) (n - p))  - log(1.0 + g * (1.0- Rsquare)) * ((double)(n - 1)));
  if (p == 1) logmargy = 0.0;
  return(logmargy);
  }

double BIC_glm(double Rsquare, int n,  int p, double SSY)
  {  double logmargy, sigma2, dn, dp;
  dn = (double) n;
  dp = (double) p;
  sigma2 = SSY*(1.0 - Rsquare);
  logmargy =  -.5*(dn*log(sigma2) +  dp*log(dn));
  return(logmargy);
  }

double AIC_glm(double Rsquare, int n,  int p, double SSY)
  {  double logmargy, sigma2, dn, dp;
  dn = (double) n;
  dp = (double) p;
  sigma2 = SSY*(1.0 - Rsquare);
  logmargy =  -.5*(dn*log(sigma2) +  dp*2.0);
  return(logmargy);
  }

double shrinkage_EB_local_glm(double R2, int n, int p, double alpha)
 {  
   /* R2 = Y'PxY/Y'Y  
      n = sample size
      p = number of rank of X (including intercept)
      alpha = prior hyperparameter
   */

   double  ghat, dn, dp, shrinkage;
   
    dn = (double) n - 1.0;
    dp = (double) p - 1.0;
    if (dp > 0.0) {
      ghat =  (((dn-dp)/dp)*R2/(1.0 - R2)) - 1.0;
      if (ghat < 0.0) ghat = 0.0;
      shrinkage = ghat/(1.0 + ghat);}
    else shrinkage = 1.0;
  return(shrinkage);
}

double logBF_hyperGprior_glm(double R2, int n, int p, double alpha)
  {  double logmargy,  a, b, c, z1, hf1;
  /* log Marginal under reference prior for phi & intercept
     hyperg prior for regression coefficients; alpha > 2 
     log marginal for null model is 0 */

  a = (double)  (n - 1) /2.0;
  b = 1.0;
  c = ((double) p - 1.0  + alpha )/2.0;
  z1 = R2;

    hf1 = hyp2f1(a, b, c, z1);
    if (p == 1) logmargy = 0.0;
    else logmargy = log(hf1) 
	     - log( (double) p - 1.0 + alpha - 2.0) + log(2.0) 
	     + log(alpha/2.0 - 1.0);
    if (! R_FINITE(logmargy))
	logmargy = logBF_hyperGprior_laplace(R2, n, p, alpha);
    return(logmargy);
  }

double logBF_hyperGprior_laplace_glm(double R2, int n, int p, double alpha)
 {  
   /* R2 = usual coefficient of determination
      n = sample size
      p = number of rank of X (including intercept)
      alpha = prior hyperparameter
      n and p are adjusted by subtrating off one 
      because of the flat prior on the intercept
   */

   double lognc, ghat, dn, dp, logmarg,sigmahat;
   
    dn = (double) n - 1.0;
    dp = (double) p - 1.0;
/*  Laplace approximation in terms of exp(tau) = g  */
    ghat = (-4.+ alpha + dp + (2. - dn)*R2 - 
	    sqrt(-8.*(-2. + alpha + dp)*(-1.0 + R2) + (-4. + alpha + dp + (2.-dn)* R2)*(-4. + alpha + dp + (2.-dn)* R2)))/(2.*(-2. + alpha + dp)*(-1. + R2)); 

    if (ghat <= 0.0)  { Rprintf("ERROR: In Laplace approximation to  logmarg,  ghat =  %f  R2 = %f p = %d  n = %d\n", ghat, R2, p,n);}
  
    sigmahat = 2.0/((dn - dp - alpha)/(( 1.0 + ghat)*(1.0 + ghat)) -
	       dn*(1.0 - R2)*(1.0 - R2)/((1.0 + ghat*(1.0 - R2))*(1.0 + ghat*(1.0 - R2))));

    /*  Substitute ghat (unconstrained) into sigma, simplify, then plug in ghat
	Seems to perform better */

    sigmahat = 2.0*(dn/((dn - dp - alpha)*(dp + alpha)))* (1.0 + ghat)*(1.0 + ghat); 
    
   sigmahat =1.0/(-ghat*(dn - alpha - dp)/(2.*(1. + ghat)*(1.+ghat)) +
     dn*(ghat*(1. - R2))/(2.*(1.+ghat*(1.-R2))*(1.+ghat*(1.-R2)))); 

    if (sigmahat <= 0 ) Rprintf("ERROR in LAPLACE APPROXIMATION to logmarg sigmhat = %f, ghat =  %f  R2 = %f p = %d  n = %d\n", sigmahat, ghat, R2, p,n); 
    lognc = log(alpha/2.0 - 1.0);
    logmarg = lognc + 
              .5*( log(2.0*PI) 
                     - (dp + alpha)*log(1.0 + ghat)
	             -  dn*log(1.0-(ghat/(1.0 + ghat))*R2)
	             + log(sigmahat)) + log(ghat);
  if (p == 1) logmarg = 0.0;
  return(logmarg);
}

double shrinkage_laplace_glm(double R2, int n, int p, double alpha)
 {  
   /* R2 = usual R2
      n = sample size
      p = number of rank of X (including intercept)
      alpha = prior hyperparameter
   */

     double  ghat, dn, dp, lognum, logden, shrinkage,sigmahat, lognc;
   
   /* numerator Laplace approximation E[g/(1 + g) | Y, model] */

    dn = (double) (n - 1);
    dp = (double) (p - 1);
    lognc = log(alpha/2.0 - 1.0);
    ghat = (-6.+ alpha + dp + (4. - dn)*R2 - 
	    sqrt(-16.*(-2. + alpha + dp)*(-1.0 + R2) + (-6. + alpha + dp + (4.-dn)* R2)*(-6. + alpha + dp + (4.-dn)* R2)))/(2.*(-2. + alpha + dp)*(-1. + R2)); 

    if (ghat <=0.0) { Rprintf("ERROR: In Laplace approximation to  E[g/(1 + g)] ghat = %f %f %d %d\n", ghat, R2, p,n ); 

           } 
   sigmahat =2.0/(ghat*(-dn + 2.+ alpha + dp)/((1. + ghat)*(1.+ghat)) +
		  dn*(ghat*(1. - R2))/((1.+ghat*(1.-R2))*(1.+ghat*(1.-R2)))); 

   if (sigmahat <= 0 ) Rprintf("ERROR in LAPLACE APPROXIMATION to E[g/(1 + g)] sigmahat = %f %f %f %d %d\n", sigmahat, ghat, R2, p,n); 

   lognum= .5*( log(2.0*PI) + 2.0*log(ghat)
		- (dp + alpha + 2.0 - dn)*log(1.0 + ghat)
		- dn*log(1.0 + ghat*(1. -R2))
		+ log(sigmahat))  +  lognc + log(ghat);
   logden = logBF_hyperGprior_laplace( R2,  n, p, alpha);
/* lognc is included here so that it cancels wth lognc for denominator
   shrinkage = exp(lognum - logden);
/*   Rprintf("%f %f %f  %f %f %f \n", ghat, sigmahat, normalpart, lognum, logden, shrinkage);  */
  if (p == 1) shrinkage = 1.0;
  return(shrinkage);
}

double log_laplace_2F1_glm(double a, double b, double c, double z)
 {  
   
   double  ghat,logmarg,sigmahat, D, A, B, C,root1, root2;

/*  Laplace approximation in terms of exp(tau) = g  */
//
//  Integrate[g^(b/2-1) (1 + g)^(a - c)/2 (1 + (1 - z) g)^(-a/2) dg 
//  Integrate[e^tau b/2 (1 + e^g)^(a - c)/2 *( 1 + (1 -z)e^g)^(-a/2) dtau  
//
   A = (b - c)*(1. - z);
   B = 2.*b - c + (a-b)*z;
   C = b;
   D = B*B - 4.*A*C;

   if (D < 0 )    Rprintf("ERROR: In Laplace approximation to 2F1");

   root1 = (-B - sqrt(D))/(2.*A);
   root2 = (-B + sqrt(D))/(2.*A);
   
   Rprintf("+root = %lf root = %lf\n", root2, root1);

   ghat = (2.*b)/(c + b*(z - 2.0) - a*z + sqrt(c*c + 4*a*b*z -  2.0*(a + b)*c *z + (a -b )*(a-b)*z*z)/(2.*(b - c)*(z - 1)));
 
   if (ghat <= 0.0)  {
     Rprintf("ERROR: In Laplace approximation to  logmarg,  ghat =  %f  z = %lf\n", ghat, z);
     ghat =  -(c +b*(-2. + z) - a*z + sqrt(c*c + 4*a*b*z -  2.0*(a + b)*c *z + (a -b )*(a-b)*z*z)/(2.*(b - c)*(z - 1.)));
   }
     sigmahat =1.0/(.5*(-a + c)*((ghat/(1 + ghat))*(1 - ghat/(1 + ghat))) +
		    a*((ghat*(1.-z))/(1.+ghat*(1.-z))*
		       (1. - (ghat*(1.-z))/(1.+ghat*(1.-z)))));

    if (sigmahat <= 0 ) Rprintf("ERROR in LAPLACE APPROXIMATION to in 2F1 sigmhat = %f, ghat =  %f  z = %f \n", sigmahat, ghat, z); 
    logmarg = .5*( log(2.0*PI) 
		   - (a - c)*log(1.0 + ghat)
		   -  a*log(1.0 + ghat*(1. - z)) 
		   + log(sigmahat)) 
               + (.5*b - 1.)*log(ghat);
  return(logmarg);
 }	


double logBF_EB_glm(double R2, int n, int p, double alpha)
 {
   double  ghat, dn, dp, logmarg;
    dn = (double) n -1.0;
    dp = (double) p -1.0;
    ghat =  (((dn-dp)/dp)*R2/(1.0 - R2)) - 1.0;
    if (ghat < 0.0) ghat = 0.0;
    logmarg = .5*( - dn*log(1.0-(ghat/(1.0 + ghat))*R2)
		   - dp*log(1.0 + ghat));
    if (p == 1) logmarg = 0.0;
    return(logmarg);
}

