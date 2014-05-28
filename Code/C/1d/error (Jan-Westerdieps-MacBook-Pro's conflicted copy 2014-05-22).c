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

int leg_Pl_e(const int l, const double x, gsl_sf_result * result) {
  //printf("%i %g %lu\n", l, x, *( long *) &x);
  if(l < 0 || x < -1.0 || x > 1.0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(l == 0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(l == 1) {
    result->val = x;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(l == 2) {
    result->val = 0.5 * (3.0*x*x - 1.0);
    result->err = GSL_DBL_EPSILON * (fabs(3.0*x*x) + 1.0);
    return GSL_SUCCESS;
  }
  else if(x == 1.0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(x == -1.0) {
    result->val = ( GSL_IS_ODD(l) ? -1.0 : 1.0 );
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(l < 100000) {
    /* upward recurrence: l P_l = (2l-1) z P_{l-1} - (l-1) P_{l-2} */

    double p_ellm2 = 1.0; /* P_0(x) */
    double p_ellm1 = x; /* P_1(x) */
    double p_ell = p_ellm1;

    double e_ellm2 = GSL_DBL_EPSILON;
    double e_ellm1 = fabs(x)*GSL_DBL_EPSILON;
    double e_ell = e_ellm1;

    int ell;

    for(ell=2; ell <= l; ell++){
      p_ell = (x*(2*ell-1)*p_ellm1 - (ell-1)*p_ellm2) / ell;
      p_ellm2 = p_ellm1;
      p_ellm1 = p_ell;

      e_ell = 0.5*(fabs(x)*(2*ell-1.0) * e_ellm1 + (ell-1.0)*e_ellm2)/ell;
      e_ellm2 = e_ellm1;
      e_ellm1 = e_ell;
    }
    result->val = p_ell;
    result->err = e_ell + l*fabs(p_ell)*GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  } else {
    GSL_ERROR ("polynomial degree too high", GSL_EDOM);
  }
}

int leg_Pl_array(const int lmax, const double x, double * result_array) {
  /* CHECK_POINTER(result_array) */

  /*
  for( int i = 0; i <= lmax; i++) {
    printf("%i %g %lu\n", i, x, *( long *) &x);
  }
  */
  if(lmax < 0 || x < -1.0 || x > 1.0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(lmax == 0) {
    result_array[0] = 1.0;
    return GSL_SUCCESS;
  }
  else if(lmax == 1) {
    result_array[0] = 1.0;
    result_array[1] = x;
    return GSL_SUCCESS;
  }
  else {
    /* upward recurrence: l P_l = (2l-1) z P_{l-1} - (l-1) P_{l-2} */

    double p_ellm2 = 1.0; /* P_0(x) */
    double p_ellm1 = x; /* P_1(x) */
    double p_ell = p_ellm1;
    int ell;

    result_array[0] = 1.0;
    result_array[1] = x;

    for(ell=2; ell <= lmax; ell++){
      p_ell = (x*(2*ell-1)*p_ellm1 - (ell-1)*p_ellm2) / ell;
      p_ellm2 = p_ellm1;
      p_ellm1 = p_ell;
      result_array[ell] = p_ell;
    }

    return GSL_SUCCESS;
  }
}

/* Given position x and degrees of freedom r, find p_r(x).
 * \[
 *   p_r(x) = \gamma_0 P_0(x) + \ldots + \gamma_{r-1} P_{r-1}(x).
 * \]
 */
double leg_sum( double *coeffs, boundary int_bounds, int r, double x) {
  double leg_res[r];
  leg_Pl_array( r-1, leg_ab2m11( int_bounds, x), leg_res);
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
  gsl_sf_result res;
  leg_Pl_e( p->n, leg_ab2m11( p->int_bounds, x), &res);
  return (p->f)(x) * res.val;
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
  double noemer = (int_bounds.b - int_bounds.a)/(2*r + 1);
  //printf("%llu gamma_%i(%llu,%i) = %g\n", l.i + r, r, l.i, l.n, res/noemer);
  return res/noemer;
}

typedef struct gamma_list {
  int i, r;
  double gamma;
  struct gamma_list *next;
} gamma_list;

gamma_list *gamma_list_create( int i, int r, double gamma) {
  gamma_list *list = malloc( sizeof( gamma_list));
  list->i = i;
  list->r = r;
  list->gamma = gamma;
  list->next = NULL;
  return list;
}

gamma_list *gamma_list_insert_after( gamma_list *list, int i, int r, double gamma) {
  gamma_list *new = gamma_list_create( i, r, gamma);
  new->next = list->next;
  list->next = new;
  return new;
}

double leg_gamma_memoize( function f, boundary int_bounds, location l, int r) {
  //return leg_gamma( f, int_bounds, l, r);
  static gamma_list *gamma[MAXN][MAXN] = {{NULL}};
  int index = (l.i + r) % MAXN; //TODO: think of a better hash function
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
  boundary int_bounds = leg_node( b, l);
  //find best polynomial
  //printf("Finding error on (%llu,%i) with degrees of freedom %i\n", l.i, l.n, r);
  //printf("%llu,%i,%i\n", l.i, l.n, r);
  assert( r > 0);
  double coeffs[r];
  int n;
  #ifdef OPENMP
  #pragma omp parallel for
  #endif
  for( n = 0; n < r; n++) {
    coeffs[n] = leg_gamma_memoize( f, int_bounds, l, n);
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
