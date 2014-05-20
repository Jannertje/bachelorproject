#include <stdio.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "macros.h"
#include "error.h"

#define ERR 1e-7

struct leg_mul_s {
  function f;
  boundary int_bounds;
  int n;
};

struct approx_quad_s {
  function f;
  double *coeffs;
  boundary int_bounds;
  int r;
};

/*
 * Linear map from [a,b] to [-1,1]. This is because Legendre Polynomials only
 * work on [-1,1].
 */
double leg_ab2m11( boundary b, double x) {
  return 2.0/(b.b-b.a) * x + (b.a+b.b)/(b.a-b.b);
}

boundary leg_node( boundary b, location l) {
  double pw2 = pow( 2.0, l.n);
  boundary res = { .a = b.a + (b.b - b.a)*l.i/pw2,
                   .b = b.a + (b.b - b.a)*(l.i+1)/pw2};
  return res;
}

/* Given position x and degrees of freedom r, find p_r(x).
 * \[
 *   p_r(x) = \gamma_0 P_0(x) + \ldots + \gamma_{r-1} P_{r-1}(x).
 * \]
 */
double leg_sum( double *coeffs, boundary int_bounds, int r, double x) {
  double leg_res[r];
  gsl_sf_legendre_Pl_array( r-1, leg_ab2m11( int_bounds, x), leg_res);
  double res = 0.0;
  int n;
  for( n = 0; n < r; n++) {
    res += coeffs[n] * leg_res[n];
  }

  return res;
}

/* For given x, find p_r(x) and
 * \[
 *   [ f(x) - p_r(x) ]^2
 * \]
 */
double leg_approx_quad( double x, void *params) {
  struct approx_quad_s *p = (struct approx_quad_s *) params;
  double leg_val = leg_sum( p->coeffs, p->int_bounds, p->r, x);
  return pow( p->f(x) - leg_val, 2.0);
}

/* 
 * For given x and n, find f(x) * P_n(x) 
 */
double leg_mul( double x, void *params) {
  struct leg_mul_s *p = (struct leg_mul_s *) params;
  return (p->f)(x) * gsl_sf_legendre_Pl( p->n, leg_ab2m11( p->int_bounds, x));
}

/*
 * \gamma_l := \int_a^b f(x) P_l(x) dx / \int_a^b P_l^2(x) dx
 * My OSX cries about gamma as a function name, so name it _gamma.
 */
double leg_gamma( function f, boundary int_bounds, location l, int n) {
  struct leg_mul_s params = { f, int_bounds, n };
  gsl_function F;
  double res, err;

  F.function = &leg_mul;
  F.params = &params;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc( 1000);
  int status = gsl_integration_qags( &F, int_bounds.a, int_bounds.b, ERR, 0.0, 1000, w, &res, &err);
  gsl_integration_workspace_free( w);
  /*
  size_t neval;
  int status = gsl_integration_qng( &F, b.a, b.b, ERR, 0.0, &res, &err, &neval);
  */
  if( status) {
    printf("hioihoi %g %i %i\n", err, status, GSL_EDIVERGE);
    exit(-1);
  }
  double noemer = (int_bounds.b - int_bounds.a)/(2*n + 1);
  //printf("gamma_%i(%i,%i) = %g\n", n, l.i, l.n, res/noemer);
  return res/noemer;
}

/* Given degrees of freedom r, find polynomial of best approximation to f on 
 * [a,b]: this is p_r(x). Then find
 * \[
 *   \| f(x) - p_r(x) \|_2^2 =: e(\node).
 * \]
 */
double error_2norm_leg( function f, boundary b, location l, int r) {
  boundary int_bounds = leg_node( b, l);
  //find best polynomial
  //printf("Finding error on (%i,%i) with degrees of freedom %i\n", l.i, l.n, r);
  printf("%i,%i,%i\n", l.i, l.n, r);
  assert( r > 0);
  double coeffs[r];
  int n;
  #ifdef OPENMP
  #pragma omp parallel for
  #endif
  for( n = 0; n < r; n++) {
    coeffs[n] = leg_gamma( f, int_bounds, l, n);
  }

  //find its two-norm
  struct approx_quad_s params = { f, coeffs, int_bounds, r };
  gsl_function F;
  double res, err;

  F.function = &leg_approx_quad;
  F.params = &params;
 
  gsl_integration_workspace *w = gsl_integration_workspace_alloc( 1000);
  int status = gsl_integration_qags( &F, int_bounds.a, int_bounds.b, ERR, 0.0, 1000, w, &res, &err);
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
