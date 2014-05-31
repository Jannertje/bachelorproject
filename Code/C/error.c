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
#define MAXG MAXN

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

typedef struct gamma_list {
  int i, r;
  double gamma;
  struct gamma_list *next;
} gamma_list;

/*
 * Linear map from [a,b] to [-1,1]. This is because Legendre Polynomials only
 * work on [-1,1].
 */
double leg_ab2m11( boundary b, double x) {
  return 2.0/(b.b-b.a) * x + (b.a+b.b)/(b.a-b.b);
}

/*
 * Construct the integration bounds from problem boundary and node location
 */
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
  double res = gsl_sf_legendre_Pl( p->n, leg_ab2m11( p->int_bounds, x));
  return (p->f)(x) * res;
}

/*
 * \gamma_r := \int_a^b f(x) P_l(x) dx / \int_a^b P_r^2(x) dx
 */
double leg_gamma( function f, boundary int_bounds, location l, int r) {
  struct leg_mul_s params = { .f = f, .int_bounds = int_bounds, .n = r };
  gsl_function F;
  double res, err;

  F.function = &leg_mul;
  F.params = &params;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc( 1000);
  int status = gsl_integration_qags( &F, int_bounds.a, int_bounds.b, ERR, 0.0, 
                                     1000, w, &res, &err);
  gsl_integration_workspace_free( w);
  if( status) {
    printf("Failed to integrate\n");
  }
  double noemer = (int_bounds.b - int_bounds.a)/(2*r + 1);
  return res/noemer;
}

gamma_list *gamma_list_create( int i, int r, double gamma) {
  gamma_list *list = malloc( sizeof( gamma_list));
  list->i = i;
  list->r = r;
  list->gamma = gamma;
  list->next = NULL;
  return list;
}

gamma_list *gamma_list_insert_after( gamma_list *list, int i, 
                                     int r, double gamma) {
  gamma_list *new = gamma_list_create( i, r, gamma);
  new->next = list->next;
  list->next = new;
  return new;
}

/*
 * Caching wrapper around leg_gamma(). Tries to locate a return value of 
 * leg_gamma for given input. If not, run leg_gamma() and put return value in 
 * hash map.
 */
double leg_gamma_memoize( function f, boundary int_bounds, location l, int r) {
  static gamma_list *gamma[MAXN][MAXG] = {{NULL}};
  int index = (l.i + r) % MAXG; //TODO: think of a better hash function?
  gamma_list *cur;
  if( (cur = gamma[l.n][index]) != NULL) {
    while( cur->next != NULL) {
      if( cur->i == l.i && cur->r == r) {
        return cur->gamma;
      }
      cur = cur->next;
    }
    if( cur->i == l.i && cur->r == r) {
      return cur->gamma;
    }

    double res = leg_gamma( f, int_bounds, l, r);
    gamma_list_insert_after( cur, l.i, r, res);
    return res;
  } else {
    double res = leg_gamma( f, int_bounds, l, r);
    gamma[l.n][index] = gamma_list_create( l.i, r, res);
    return res;
  }
}

/* Given degrees of freedom r, find polynomial of best approximation to f on 
 * [a,b]: this is p_r(x). Then find
 * \[
 *   \| f(x) - p_r(x) \|_2^2 =: e(\node).
 * \]
 */
double error_2norm_leg( function f, boundary b, location l, int r) {
  //construct int_bounds once, possibly costly operation
  boundary int_bounds = leg_node( b, l);
  assert( r > 0);
  double coeffs[r];
  int n;

  //construct polynomial coefficients
  for( n = 0; n < r; n++) {
    coeffs[n] = leg_gamma_memoize( f, int_bounds, l, n);
  }

  //find the 2-norm of the error
  struct approx_quad_s params = { f, coeffs, int_bounds, r };
  gsl_function F;
  double res, err;

  F.function = &leg_approx_quad;
  F.params = &params;
 
  gsl_integration_workspace *w = gsl_integration_workspace_alloc( 1000);
  int status = gsl_integration_qags( &F, int_bounds.a, int_bounds.b, ERR, 0.0, 
                                     1000, w, &res, &err);
  gsl_integration_workspace_free( w);
  if( status) {
    printf("Failed to integrate\n");
    exit(-1);
  }
  return res;
}
