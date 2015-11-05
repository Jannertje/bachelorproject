#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_errno.h>
#include <assert.h>
#include <stdio.h>
#include "types.h"
#include "workspace.h"
#include "tree.h"
#include "helper.h"
#include "main.h"

/* GSL struct which is passed to gamma_x */
typedef struct gamma_x_p {
  double a, b;
  int j;
} gamma_x_p;

/* simple wrapper around gsl_integrate_* */
double integrate( double (* function)( double x, void *params), void *params,
                  double a, double b, double epsabs, double epsrel) {
  double res, err;
  gsl_set_error_handler_off();
  gsl_function F = {.function = function, .params = params};
  gsl_integration_workspace *w = gsl_integration_workspace_alloc( 10000);
  gsl_integration_qags( &F, a, b, epsabs, epsrel, 10000, w, &res, &err);
  gsl_integration_workspace_free( w);
  /*
size_t neval;
  gsl_integration_qng( &F, a, b, ERR, 0.0, &res, &err, &neval);
  */
  return res;
}

double gamma_x( double x, void *params) {
  gamma_x_p *p = (gamma_x_p *) params;
  return f(x) * gsl_sf_legendre_Pl( p->j, 2.0/(p->b-p->a)*x+(p->a+p->b)/(p->a-p->b));
}

double norm_f( double x, void *params) {
  double fx = f(x);
  return fx*fx;
}

double construct_p_n( tree *node, int n) {
  int j;
  double gamma, sum_gamma_sq = 0;
  //make sure there is space in the gamma array
  tree_assert_gamma_len( node, n);
  #pragma omp parallel for private(j, gamma) reduction(+:sum_gamma_sq)
  for( j = 0; j < n; j++) {
    gamma = tree_get_gamma( node, j);
    if( gamma != gamma) { //is gamma set?
      gamma_x_p params = { .a = node->a, .b = node->b, .j = j};
      gamma = (2.0*j+1.0)/(node->b - node->a) * 
              integrate( &gamma_x, &params, node->a, node->b, 1E-15, 0.0);
      #pragma omp critical
      {
        tree_set_gamma( node, j, gamma);
      }
    }
    double norm_phi_j = (node->b - node->a)/(2.0*j+1.0);
    sum_gamma_sq += gamma*gamma*norm_phi_j;
    //fprintf( stderr, "gamma_%i([%g,%g]) = %g\n", j, node->a, node->b, gamma);
  }
  return sum_gamma_sq;
}

double error( tree *node, int r) {
  assert( r > 0);
  double norm_p_n_sq = construct_p_n( node, r);

  double e = integrate( &norm_f, NULL, node->a, node->b, 1E-15, 0.0);
  //fprintf( stderr, "error_%i([%g,%g]) = %g - %g = %g\n", r, node->a, node->b, e, norm_p_n_sq, e - norm_p_n_sq);
  return fmax( 0, e - norm_p_n_sq);
}
