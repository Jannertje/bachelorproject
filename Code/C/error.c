#include <gsl/gsl_integration.h>
#include <stdio.h>
#include "tri.h"
#include "pkd.h"
#include "pair.h"

#define ERR 1E-5

typedef struct gamma_s {
  tri *node;
  double norm;
  int j, k;
} gamma_s;

typedef struct gamma_y_s {
  double x;
  gamma_s *s;
} gamma_y_s;

typedef struct error_s {
  tri *node;
  double *coeffs;
  int l;
} error_s;

typedef struct error_y_s {
  error_s *s;
  double x;
} error_y_s;

point a = {0, 0};
point b = {1, 0};
point c = {0, 1};
tri *node;

double f( double x, double y) {
  return 1;
}

/* simple wrapper around gsl_integrate_* */
double integrate( double (* function)( double x, void *params), void *params,
                  double a, double b) {
  gsl_function F = {.function = function, .params = params};
  /*
  gsl_integration_workspace *w = gsl_integration_workspace_alloc( 1000);
  gsl_integration_qags( &F, -1.0, x, ERR, 0.0, 1000, w, &res, &err);
  gsl_integration_workspace_free( w);
  */
  double res, err;
  size_t neval;
  gsl_integration_qng( &F, a, b, ERR, 0.0, &res, &err, &neval);
  return res;
}

/* compute f(G(x,y)) * Q_{j,k}(x,y),
 * in other words, the integrand of gamma_{j,k}. */
double gamma_y( double y, void *params) {
  gamma_y_s *p = (gamma_y_s *) params;
  point q = {.x = p->x, .y = y};
  point r = tri_ref2node( p->s->node, q);
  point Fr = {.x = r.x*(1-r.y), r.y};
  point Fq = {.x = q.x*(1-q.y), q.y};
  double fr = f(r.x, r.y);
  double Qp = pkd_eval( p->s->j, p->s->k, q);
  return fr * Qp*(1-y);
}

/* compute \int_0^{1-x} gamma_y( x, y) dy,
 * in other words, the inner integral of gamma_{j,k}. */
double gamma_x( double x, void *params) {
  gamma_s *p = (gamma_s *) params;
  gamma_y_s p2 = {.s = p, .x = x};
  return integrate( &gamma_y, &p2, 0, 1);
}

/* compute \int_0^1 gamma_x( x, y) dx,
 * in other words, find the outer integral of gamma_{j,k}
 * which is needed for the polynomial of best approximation */
double error_gamma( int j, int k) {
  double norm = pkd_norm( j, k);
  gamma_s params = { .node = node, .norm = norm, .j = j, .k = k};
  double n = pkd_norm( j, k);
  return 2*node->vol*integrate( &gamma_x, &params, 0, 1)/(n*n);
}

//TODO: extend to other triangles
/* compute [f(x,y) - p_n(G^{-1}(x,y))]^2 */
double error_y( double y, void *params) {
  error_y_s *p = (error_y_s *) params;
  point q = {.x = p->x, .y = y};
  point r = tri_ref2node( p->s->node, q);

  double sum = 0.0;
  int j, k, l;
  for( l = 0; l < p->s->l; l++) {
    pair_invcantor( &j, &k, l);
    sum += p->s->coeffs[l] * pkd_eval( j, k, r);
  }

  double err = f(q.x, q.y) - sum;
  return err*err;
}

/* compute \int_0^{1-x} error_y(x, y) dy */
double error_x( double x, void *params) {
  //integrate
  error_s *p = (error_s *) params;
  error_y_s p2 = {.s = p, .x = x};
  return integrate( &error_y, &p2, 0, 1.0-x);
}

/* find gamma_{j,k} for j \in \{0, \ldots, n\} and k \in \{0, \ldots, n-j\} */
void construct_p_n( int n, double *coeffs) {
  int j, k, l;
  for( j = 0; j < n; j++) {
    for( k = 0; k < n-j; k++) {
      pair_cantor( j, k, &l);
      coeffs[l] = error_gamma( j, k);
      printf("%i %i %g, %g\n", j, k, coeffs[l], pkd_norm( j, k));
    }
  }
}

/* compute e_n( \node) by computing 
 * \|f - p_n\|^2_{2,\node} = \int_0^1 error_x( x, y) dx. */
double error( int n) {
  //find polynomial
  double coeffs[(n+2)*(n+1)/2];
  construct_p_n( n, coeffs);

  //integrate
  error_s params = { .node = node, .coeffs = coeffs, .l = (n+2)*(n+1)/2};
  return integrate( &error_x, &params, 0, 1);
}

/*
int main() {
  node = tri_create( a, b, c);
  double g = error( 5);

  //double g = error_gamma( 1, 1);
  printf("%g\n", g);
  return 0;
}
*/
