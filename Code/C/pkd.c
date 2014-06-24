#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tri.h"

double pkd_norm( int j, int k) {
  return sqrt(2.0/((2*j+1)*(j+k+1)));
}

//find P_k^(a, 0)(x) for x \in [-1,1]
//method same as legendre
double jacobi( int k, double a, double x) {
  if(k < 0 || x < -1.0 || x > 1.0) {
    printf("domain error: %i, %g\n", k, x);
    exit(1);
  }
  else if(k == 0) {
    return 1.0;
  }
  else if(k == 1) {
    return ((a+2)*x+a)/2;
  }
  else if(k < 100000) {
    //hier zit de fout
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
  else {
    printf("Order too high!\n");
    exit(1);
  }
}

double pkd_eval_square( int j, int k, point p) {
  double leg = gsl_sf_legendre_Pl( j, p.x);
  double pwr = gsl_pow_int( (1.0-p.y)/2.0, j);
  double jac = jacobi( k, 2.0*j+1.0, p.y);
  return leg*pwr*jac;
}

/*
int main() {
  point a = {-1, -1};
  point b = {1, -1};
  point c = {-1, 1};
  tri *node = tri_create( a, b, c);
  point p = {0.33,0.75};
  double Qp = pkd_eval( node, 5, 5, p);
  return 0;
}
*/
