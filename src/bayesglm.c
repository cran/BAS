#include "sampling.h"
#include "family.h"
#include <R_ext/BLAS.h>

typedef struct glmfamilystruc {
  const char *family;
  const char *link;
  void (*mu_eta)(double *eta, double *mu, int n);
  void (*linkfun)(double *mu, double *eta, int n);
  void (*variance)(double * mu, double *var, int n);
  void (*dev_resids)(double *y, double *mu, double *weights, double *resids, int n);
  void (*linkinv)(double *eta, double *mu, int n);
  double (*dispersion)(double *resid,  double *weights, int n, int rank);
} glmstptr;

typedef struct coefpriorstruc {
  const char *family;
  const char *class;
  double *hyper;
  double (*log_marginal_likelihood)(double dev, double regSS, int n, int p, int pgamma, double g, double *hyper);
  double (*shrinkage)(double dev,  double regSS, int n, int p, int pgamma, double g, double *hyper);
  double (*g)(double dev,  double regSS, int n, int p, int pgamma, double *hyper);
} coefdistptr;




 double no_shrinkage(double dev,  double regSS, int n, int p, int pgamma, double g,  double *hyper) {
  return 1.0;
}

 double shrinkage_gprior(double dev,  double regSS, int n, int p, int pgamma,  double g, double *hyper) {
  return g/(1.0 + g) ;
}

 double g_EB_local(double dev,  double regSS, int n, int p, int pgamma, double *hyper) { 
  double g = regSS/pgamma - 1.0;
  return g > 0.0 ? g : 0.0;
}

 double g_gprior(double dev,  double regSS, int n, int p, int pgamma, double *hyper) { 
  return hyper[0];
}

 double no_g(double dev,  double regSS, int n, int p, int pgamma, double *hyper) { 
  return 1.0;
}

double log_marginal_likelihood_IC(double dev, double regSS, int n, int p, int pgamma, double g, double *hyper) {
  return -.5*(dev +  pgamma*hyper[0]);
}
double log_marginal_likelihood_gprior(double dev, double regSS, int n, int p, int pgamma, double g, double *hyper) {
  return -.5*(dev  + pgamma*log(g + 1.0) + regSS/(g + 1.0));
}


/* Version of glm.fit that can be called directly from R or C*/

SEXP glm_fit(SEXP RX, SEXP RY,SEXP family, SEXP Roffset, SEXP Rweights, SEXP Rstarteta, SEXP Rstartcoef, SEXP Rpriorcoef, SEXP Repsilon)
{
  int   *xdims = INTEGER(getAttrib(RX,R_DimSymbol)), n=xdims[0], p = xdims[1];
  int inc = 1,  job = 01, info = 1, nmodels=1;
  
  SEXP ANS = PROTECT(allocVector(VECSXP, 10));
  SEXP RXwork = PROTECT(duplicate(RX)), 
    RYwork = PROTECT(duplicate(RY)),  
    RWwork = PROTECT(duplicate(RY)), Rvariance = PROTECT(duplicate(RY)), 
    Rmu_eta = PROTECT(duplicate(RY)), Reta= PROTECT(duplicate(Rstarteta)),  
    Rmu = PROTECT(duplicate(RY)),
    Rcoef= PROTECT(duplicate(Rstartcoef)),
    Rcoefwork = PROTECT(duplicate(Rstartcoef)), Rrank=PROTECT(allocVector(INTSXP,1)),
    Rresdf = Rrank=PROTECT(allocVector(INTSXP,1)),
    Rcov = PROTECT(allocVector(REALSXP, p*p)),      RR = PROTECT(allocVector(REALSXP, p*p)),  
    Rse= PROTECT(allocVector(REALSXP, p)),  
    Rresiduals= PROTECT(duplicate(RY)), Reffects= PROTECT(duplicate(RY)),
    Rpivot=PROTECT(allocVector(INTSXP,p)),
    Rqrauxmat=PROTECT(allocVector(REALSXP,p)), 
    Rworkmat=PROTECT(allocVector(REALSXP,2*p)), Rlog_marg_lik=PROTECT(allocVector(REALSXP,nmodels)), 
    Rdeviance=PROTECT(allocVector(REALSXP,nmodels)), RregSS=PROTECT(allocVector(REALSXP,nmodels)), 
    Rg = PROTECT(allocVector(REALSXP, nmodels)), Rshrinkage = PROTECT(allocVector(REALSXP, nmodels));


  double *X=REAL(RX), *Y=REAL(RY), *Xwork=REAL(RXwork),
    *w=REAL(RWwork),*Ywork=REAL(RYwork), *effects=REAL(Reffects),
    *coef=REAL(Rcoef),*coefwork=REAL(Rcoefwork), *se=REAL(Rse), *cov = REAL(Rcov), *R = REAL(RR),
    *work=REAL(Rworkmat), *qraux=REAL(Rqrauxmat), *weights=REAL(Rweights),
    *mu=REAL(Rmu), *offset=REAL(Roffset),*eta=REAL(Reta),  *mu_eta=REAL(Rmu_eta),
    *residuals=REAL(Rresiduals), *dev=REAL(Rdeviance), *regSS = REAL(RregSS), *g = REAL(Rg),
    *variance=REAL(Rvariance), *hyper;

  double  zero = 0.0, one = 1.0,  tol = DBL_EPSILON, devold, devnew, regss;

  int   i, j, l, m, rank=1, *pivot=INTEGER(Rpivot), conv=0;

  glmstptr *glmfamily;
  coefdistptr *coefprior;
  char uplo[] = "U", trans[]="N";
  
  coefprior = (struct coefpriorstruc *) R_alloc(1, sizeof(struct coefpriorstruc));
  coefprior->family = CHAR(STRING_ELT(getListElement(Rpriorcoef, "family"),0));
  coefprior->class  = CHAR(STRING_ELT(getListElement(Rpriorcoef, "class"),0));
  if (getListElement(Rpriorcoef, "hyper") != R_NilValue) 	
    hyper = REAL(getListElement(Rpriorcoef, "hyper"));	

  if  (strcmp(coefprior->class, "gprior") == 0) {
    coefprior->shrinkage = shrinkage_gprior;
    coefprior->log_marginal_likelihood = log_marginal_likelihood_gprior;
    if  (strcmp(coefprior->family, "fixed-g-prior") == 0) 
      coefprior->g = g_gprior;
    else 
      coefprior->g = g_EB_local;
}
  if  (strcmp(coefprior->class, "IC") == 0) {
    coefprior->shrinkage = no_shrinkage;
    coefprior->log_marginal_likelihood = log_marginal_likelihood_IC;
    coefprior->g = no_g;
}
  
  Rprintf("prior %s\n", coefprior->family);

  glmfamily = (struct glmfamilystruc *) R_alloc(1, sizeof(struct glmfamilystruc));
  glmfamily->family = CHAR(STRING_ELT(getListElement(family, "family"),0));
  glmfamily->link = CHAR(STRING_ELT(getListElement(family, "link"),0));
  
  Rprintf("link %s\n", glmfamily->link);
  if  (strcmp(glmfamily->family, "binomial") == 0) {
    glmfamily->dev_resids = binomial_dev_resids;
    glmfamily->dispersion = binomial_dispersion;
    if (strcmp(glmfamily->link, "logit") == 0) {
       glmfamily->linkfun = logit_link;	
       glmfamily->mu_eta = logit_mu_eta;
       glmfamily->variance = logit_variance; 
       glmfamily->linkinv =  logit_linkinv;
	}
   else  Rprintf("no other links implemented yet\n");
  }
  else  Rprintf("no other families implemented yet\n");
   
   

  for (m=0; m< nmodels; m++){
    glmfamily->linkinv(eta, mu, n);
    glmfamily->dev_resids(Y, mu, weights, residuals, n);
    devold = deviance(residuals, n);
    devnew = devold;

  while ( conv < 1) {
    for (j=0; j<p; j++) {
      pivot[j] = j+1;
  }	

    glmfamily->mu_eta(eta, mu_eta, n);
    glmfamily->variance(mu, variance, n);

    for (i=0, l=0; i<n; i++) {
      w[i] = sqrt(weights[i]*mu_eta[i]*mu_eta[i]/variance[i]);
      Ywork[i] = w[i]*(eta[i] - offset[i] + (Y[i] - mu[i])/mu_eta[i]);
      residuals[i] = (Y[i] - mu[i])/mu_eta[i];
    }
    for (j=0, l=0; j<p; j++) {
      for (i=0; i<n; i++, l++) {
	Xwork[l] = REAL(RX)[l]*w[i];
      }
    }

    rank = 1;

    F77_NAME(dqrls)(&Xwork[0], &n, &p, &Ywork[0], &inc, &tol,  &coefwork[0],
	 &residuals[0], &effects[0], &rank, &pivot[0], &qraux[0], &work[0]);
    
    for (j=0; j<p; j++) { 
      coef[pivot[j] - 1] = coefwork[j];}

    F77_NAME(dcopy)(&n, &offset[0], &inc, &eta[0], &inc);
    F77_NAME(dgemv)(trans, &n, &p, &one, &X[0], &n, &coef[0], &inc, &one, &eta[0],&inc);

    glmfamily->linkinv(eta, mu, n);
    glmfamily->dev_resids(Y, mu, weights, residuals, n);
    devnew = deviance(residuals, n);
    glmfamily->mu_eta(eta, mu_eta, n);
    glmfamily->variance(mu, variance, n);

    devnew = deviance(residuals, n);
    //    Rprintf("old %f new %f conv %d\n", devold,devnew, conv);

    if (fabs(devnew - devold)/(0.1 + fabs(devnew)) < REAL(Repsilon)[0]) {
     conv = 1;
	  }
    else { devold=devnew;}
  }

  dev[m] = devnew;

  chol2se(&Xwork[0], &se[0], &R[0], &cov[0], p, n);
  Rprintf("compute marginal lik quantitis");
  regSS[m] = quadform(coef, coefwork, R, p);
  g[m] = coefprior->g(dev[m],  regSS[m],  n,  p,  rank, hyper);
  REAL(Rshrinkage)[m]  = coefprior->shrinkage(dev[m],  regSS[m],  n,  p, rank, g[m],  hyper);
  REAL(Rlog_marg_lik)[m]  = coefprior->log_marginal_likelihood(dev[m],  regSS[m],  n,  p, rank, g[m],  hyper);

  INTEGER(Rrank)[m] = rank;
  }

  SET_VECTOR_ELT(ANS, 0, Rcoef);
  SET_VECTOR_ELT(ANS, 1, Rse);
  SET_VECTOR_ELT(ANS, 2, Rmu);
  SET_VECTOR_ELT(ANS, 3, Rdeviance);
  SET_VECTOR_ELT(ANS, 4, Rrank);
  SET_VECTOR_ELT(ANS, 5, Rg);
  SET_VECTOR_ELT(ANS, 6, Rshrinkage);
  SET_VECTOR_ELT(ANS, 7, RregSS);
  SET_VECTOR_ELT(ANS, 8, Rlog_marg_lik);

  // SET_VECTOR_ELT(ANS, 5, Rresdf);
  
  UNPROTECT(25);
  
  return(ANS);
  }


