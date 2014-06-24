#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <assert.h>
#include <stdio.h>
#include "tri.h"
#include "types.h"
#include "workspace.h"
#include "partition.h"
#include "tree.h"
#include "pkd.h"
#include "pair.h"
#include "main.h"

#define ERR 1E-5

typedef struct gamma_s_p {
  workspace *w;
  tri *t;
  double norm;
  int j, k;
} gamma_s_p;

typedef struct gamma_t_p {
  double s;
  gamma_s_p *p;
} gamma_t_p;

typedef struct error_s_p {
  workspace *w;
  tri *t;
  double *coeffs;
  int l;
} error_s_p;

typedef struct error_t_p {
  double s;
  error_s_p *p;
} error_t_p;

point F( point st) {
  point Fst = {.x = (1.0-st.y)*(1.0+st.x)/4.0, .y = (st.y + 1.0)/2.0};
  return Fst;
}

/* simple wrapper around gsl_integrate_* */
double integrate( double (* function)( double x, void *params), void *params,
                  double a, double b) {
  double res, err;
  gsl_set_error_handler_off();
  gsl_function F = {.function = function, .params = params};
  gsl_integration_workspace *w = gsl_integration_workspace_alloc( 1000);
  gsl_integration_qags( &F, a, b, ERR, 0.0, 1000, w, &res, &err);
  gsl_integration_workspace_free( w);
  /*
  size_t neval;
  gsl_integration_qng( &F, a, b, ERR, 0.0, &res, &err, &neval);
  */
  return res;
}

/* compute (1-t)*f(G(F(s,t)))*Q_{j,k}^\square(s,t)
 * in other words, the integrand of gamma_{j,k}. */
double gamma_t( double t, void *params) {
  gamma_t_p *p = (gamma_t_p *) params;
  point st = {.x = p->s, .y = t};
  point Fst = F( st);
  point GFst = tri_ref2t( p->p->w, p->p->t, Fst);
  double fr = f(GFst.x, GFst.y);
  double Qp = pkd_eval_square( p->p->j, p->p->k, st);
  return fr * Qp*(1-t);
}

/* compute \int_{-1}^1 gamma_y( s, t) dt,
 * in other words, the inner integral of gamma_{j,k}. */
double gamma_s( double s, void *params) {
  gamma_s_p *p = (gamma_s_p *) params;
  gamma_t_p p2 = {.p = p, .s = s};
  return integrate( &gamma_t, &p2, -1, 1);
}

/* compute \int_0^1 gamma_x( x, y) dx,
 * in other words, find the outer integral of gamma_{j,k}
 * which is needed for the polynomial of best approximation */
double error_gamma( workspace *w, tri *t, int j, int k) {
  double norm = pkd_norm( j, k);
  gamma_s_p params = { .w = w, .t = t, .norm = norm, .j = j, .k = k};
  fprintf( stderr, "gamma %i %i %i: %i %i\n", t->p[0], t->p[1], t->p[2], j, k);
  return (2.0*j+1.0)*(j+k+1.0)/4.0 * integrate( &gamma_s, &params, -1, 1);
}

/* compute (1-t)[f(G(F(s,t))) - \sum \sum \gamma_{j,k} Q_{j,k}^\square(s,t)]^2 */
double error_t( double t, void *params) {
  error_t_p *p = (error_t_p *) params;
  point st = {.x = p->s, .y = t};
  point Fst = F( st);
  point GFst = tri_ref2t( p->p->w, p->p->t, Fst);

  double sum = 0.0;
  int j, k, l;
  for( l = 0; l < p->p->l; l++) {
    pair_invcantor( &j, &k, l);
    sum += p->p->coeffs[l] * pkd_eval_square( j, k, st);
  }

  double err = f(GFst.x, GFst.y) - sum;
  return (1-t)*err*err;
}

/* compute \int_0^{1-x} error_y(x, y) dy */
double error_s( double s, void *params) {
  error_s_p *p = (error_s_p *) params;
  error_t_p p2 = {.p = p, .s = s};
  return integrate( &error_t, &p2, -1, 1);
}

/* find gamma_{j,k} for j \in \{0, \ldots, n\} and k \in \{0, \ldots, n-j\} */
void construct_p_n( workspace *w, tri *t, int n, double *coeffs) {
  int j, k, l;
  for( j = 0; j < n; j++) {
    for( k = 0; k < n-j; k++) {
      pair_cantor( j, k, &l);
      coeffs[l] = error_gamma( w, t, j, k);
      //printf("gamma: %i %i %g\n", j, k, coeffs[l]);
    }
  }
}

/* compute e_n( \node) by computing 
 * \|f - p_n\|^2_{2,\node}
 */
double error( workspace *w, tree *node, int n) {
  assert( n > 0);
  //find polynomial
  double *coeffs = malloc((n+1)*n/2 * sizeof( double));
  construct_p_n( w, w->tris[node->i], n, coeffs);

  if( node->hp) {
    hptree_set_coeffs( w, node, n, coeffs);
  } else {
    node->info.h->coeffs = coeffs;
  }

  //integrate
  tri *t = w->tris[node->i];
  error_s_p params = { .w = w, .t = w->tris[node->i], 
                       .coeffs = coeffs, .l = (n+1)*n/2};
  double e = t->vol/4.0*integrate( &error_s, &params, -1, 1);
  //fprintf( stderr, "error %i\n", n);
  return e;
}

/*
int main( void) {
  int i;
  workspace *w = workspace_init();
  point points[4] = {
    {.x = 0, .y = 0}, 
    {.x = 1, .y = 0}, 
    {.x = 0, .y = 1}, 
    {.x = 1, .y = 1}
  };

  for( i = 0; i < 4; i++) {
    workspace_add_point( w, points[i]);
  }

  tri *tris[2] = {
    tri_create( w, 0, 1, 2), 
    tri_create( w, 3, 2, 1)
  };

  for( i = 0; i < 2; i++) {
    workspace_add_tri( w, tris[i]);
  }

  partition_setup( w);

  double g = error( w, w->tris[w->leaves[0]->i], 5);

  //double g = error_gamma( 1, 1);
  printf("%g\n", g);
  return 0;
}
*/
