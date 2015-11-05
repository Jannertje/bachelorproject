#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tri.h"

/* See Definition 2.14. Find P_k^(a,0)(x) by upward recurrence of 
 * Theorem 2.15 */
double jacobi( int k, double a, double x) {
  if(k < 0 || x < -1.0 || x > 1.0) {
    printf("domain error: %i, %g\n", k, x);
    exit(1);
  } else if(k == 0) {
    return 1.0;
  } else if(k == 1) {
    return ((a+2)*x+a)/2;
  } else {
    /* upward recurrence */

    double p_jm2 = 1.0;
    double p_jm1 = ((a+2)*x+a)/2;
    double p_j = p_jm1;

    int j;

    for(j=2; j <= k; j++){
      p_j = (2.0*j+a-1.0)*((2.0*j+a)*(2.0*j+a-2.0)*x + a*a)*p_jm1 
          - 2.0*(j+a-1.0)*(j-1.0)*(2.0*j+a)*p_jm2;
      p_j = p_j / (2.0*j*(j+a)*(2.0*j+a-2.0));
      p_jm2 = p_jm1;
      p_jm1 = p_j;
    }

    return p_j;
  }
}

/* See Remark 2.17 */
double pkd_eval_square( int j, int k, point p) {
  double leg = gsl_sf_legendre_Pl( j, p.x);
  double pwr = gsl_pow_int( (1.0-p.y)/2.0, j);
  double jac = jacobi( k, 2.0*j+1.0, p.y);
  return leg*pwr*jac;
}
