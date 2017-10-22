#include <R.h>
#include <Rmath.h>





/******* Sampling truncated Gamma distribution *********
 alpha, beta: parameters of Gamma distribution
 a, b: lower and upper bound
 n: number of samples to generate

 reference: Damien and Walker 2001
 */

void r_trunc_gamma(int *n, double *alpha, double *beta, double *a, double *b, double *x){
  int i;
  // w = beta * x, aw and bw are lower and upper bound of w
  double y, z, aw = (*a) * (*beta), bw = (*b) * (*beta), w = (*a + *b) / 2 * (*beta);

  for(i=0; i< *n; i++){
    y = runif(0, exp(-w));
    z = runif(0, 1);
    if(bw < -log(y)){
      w = R_pow((z * R_pow(bw, *alpha) + (1 - z) * R_pow(aw, *alpha)), 1 / *alpha);
    }
    else {
      w = R_pow((z * R_pow(-log(y), *alpha) + (1 - z) * R_pow(aw, *alpha)), 1 / *alpha);
    }

    x[i] = w / (*beta);
  }

}


