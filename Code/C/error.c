#include <stdio.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "error.h"

#define ERR 1e-7

struct leg_mul_s {
  function f;
  boundary b;
  int l;
};

struct approx_quad_s {
  function f;
  double *coeffs;
  boundary b;
  int r;
};

/*
 * Linear map from [a,b] to [-1,1]. This is because Legendre Polynomials only
 * work on [-1,1].
 */
double ab2m11( boundary b, double x) {
  double res = 2.0/(b.b-b.a) * x + (b.a+b.b)/(b.a-b.b);
  return res;
}

/* 
 * For given x and l, find f(x) * P_l(x) 
 */
double leg_mul( double x, void *params) {
  struct leg_mul_s *p = (struct leg_mul_s *) params;
  return (p->f)(x) * gsl_sf_legendre_Pl( p->l, ab2m11( p->b, x));
}

/*
 * \gamma_l := \int_a^b f(x) P_l(x) dx / \int_a^b P_l^2(x) dx
 * My OSX cries about gamma as a function name, so name it _gamma.
 */
double _gamma( function f, boundary b, int l) {
  struct leg_mul_s params = { f, b, l };
  gsl_function F;
  double teller, err;

  F.function = &leg_mul;
  F.params = &params;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc( 1000);
  int status = gsl_integration_qags( &F, b.a, b.b, ERR, 0.0, 1000, w, &teller, &err);
  gsl_integration_workspace_free( w);
  /*
  size_t neval;
  int status = gsl_integration_qng( &F, b.a, b.b, ERR, 0.0, &teller, &err, &neval);
  */
  if( status) {
    printf("hioihoi %g %i %i\n", err, status, GSL_EDIVERGE);
    exit(-1);
  }
  double noemer = (b.b-b.a)/(2*l + 1);
  printf("gamma_%i[%f,%f] = %g\n", l, b.a, b.b, teller/noemer);
  return teller/noemer;
}

/* Given position x and degrees of freedom r, find p_r(x).
 * \[
 *   p_r(x) = \gamma_0 P_0(x) + \ldots + \gamma_{r-1} P_{r-1}(x).
 * \]
 */
double leg_sum( double *coeffs, boundary b, int r, double x) {
  double leg_res[r];
  gsl_sf_legendre_Pl_array( r-1, ab2m11( b, x), leg_res);
  double res = 0.0;
  int l;
  for( l = 0; l < r; l++) {
    res += coeffs[l] * leg_res[l];
  }

  return res;
}

/* For given x, find
 * \[
 *   [ f(x) - p_r(x) ]^2
 * \]
 */
double approx_quad( double x, void *params) {
  struct approx_quad_s *p = (struct approx_quad_s *) params;
  double leg_val = leg_sum( p->coeffs, p->b, p->r, x);
  return pow( p->f(x) - leg_val, 2.0);
}

/* Given degrees of freedom r, find polynomial of best approximation to f on 
 * [a,b]: this is p_r(x). Then find
 * \[
 *   \| f(x) - p_r(x) \|_2^2 =: e(\node).
 * \]
 */
double error_2norm( function f, boundary b, int r) {
  //find best polynomial
  //printf("Finding error on [%g,%g] with degrees of freedom %i\n", b.a, b.b, r);
  assert( r > 0);
  double coeffs[r];
  int l;
  #ifdef OPENMP
  #pragma omp parallel for
  #endif
  for( l = 0; l < r; l++) {
    coeffs[l] = _gamma( f, b, l);
  }

  //find its two-norm
  struct approx_quad_s params = { f, coeffs, b, r };
  gsl_function F;
  double res, err;

  F.function = &approx_quad;
  F.params = &params;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc( 1000);
  int status = gsl_integration_qags( &F, b.a, b.b, ERR, 0.0, 1000, w, &res, &err);
  gsl_integration_workspace_free( w);
  /*
  size_t neval;
  int status = gsl_integration_qng( &F, b.a, b.b, ERR, 0.0, &res, &err, &neval);
  */
  if( status) {
    printf("hioihoi2\n");
    exit(-1);
  }
  return res;
}
