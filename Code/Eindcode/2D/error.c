#include <omp.h>
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
#include "helper.h"
#include "pair.h"
#include "main.h"

#define ERR 1E-5
int omp_thread_count() {
    int n = 0;
    #pragma omp parallel reduction(+:n)
    n += 1;
    return n;
}

/* GSL struct which is passed to gamma_s */
typedef struct gamma_s_p {
  workspace *w;
  tri *t;
  int j, k;
} gamma_s_p;

/* GSL struct which is passed to gamma_t */
typedef struct gamma_t_p {
  double s;
  gamma_s_p *p;
} gamma_t_p;

/* GSL struct which is passed to norm_f_s */
typedef struct norm_f_s_p {
  workspace *w;
  tri *t;
} norm_f_s_p;

/* GSL struct which is passed to norm_f_t */
typedef struct norm_f_t_p {
  double s;
  norm_f_s_p *p;
} norm_f_t_p;

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

/* compute \int_{-1}^1 gamma_s( s, t) ds,
 * in other words, find the outer integral of gamma_{j,k}
 * which is needed for the polynomial of best approximation */
double error_gamma( workspace *w, tri *t, int j, int k) {
  gamma_s_p params = { .w = w, .t = t, .j = j, .k = k};
  //fprintf( stderr, "gamma %i %i %i: %i %i\n", t->p[0], t->p[1], t->p[2], j, k);
  return (2.0*j+1.0)*(j+k+1.0)/4.0 * integrate( &gamma_s, &params, -1, 1);
}

/* compute (1-t)[f(G(F(s,t)))]^2 to find the norm of f */
double norm_f_t( double t, void *params) {
  norm_f_t_p *p = (norm_f_t_p *) params;
  point st = {.x = p->s, .y = t};
  point Fst = F( st);
  point GFst = tri_ref2t( p->p->w, p->p->t, Fst);

  double err = f(GFst.x, GFst.y);
  return (1-t)*err*err;
}

/* compute \int_{-1}^1 norm_f_t(s, t) dt to find the norm of f */
double norm_f_s( double s, void *params) {
  norm_f_s_p *p = (norm_f_s_p *) params;
  norm_f_t_p p2 = {.p = p, .s = s};
  return integrate( &norm_f_t, &p2, -1, 1);
}

/* find gamma_{j,k} for j \in \{0, \ldots, n\} and k \in \{0, \ldots, n-j\}.
 * at the same time, construct the value of Theorem 2.29 */
double construct_p_n( workspace *w, tree *node, int n) {
  int j, k, l;
  double gamma, sum_gamma_sq = 0;
  //make sure there is space in the gamma array
  tree_assert_gamma_len( w, node, n);
  #pragma omp parallel for private(k, l, gamma)
  for( j = node->hgammas; j < n; j++) {
    pair_invcantor( &k, &l, j);
    #if CACHE_GAMMA
      gamma = tree_get_gamma( w, node, j);
      if( gamma != gamma) { //is gamma set?
        gamma = error_gamma( w, w->tris[node->i], k, l);
        #pragma omp critical
        {
          tree_set_gamma( w, node, j, gamma);
        }
      }
    #else
      gamma = error_gamma( w, w->tris[node->i], k, l);
      #pragma omp critical
      {
        tree_set_gamma( w, node, j, gamma);
      }
    #endif
    double norm_phi_k_l = w->tris[node->i]->vol / ((2.0*k + 1.0)*(k+l+1.0));
    sum_gamma_sq += gamma*gamma*norm_phi_k_l;
  }
  return sum_gamma_sq;
}

double find_norm_f( workspace *w, tree *node) {
  #if CACHE_NORM_F
    if( node->norm_f == -1) {
      tri *t = w->tris[node->i];
      norm_f_s_p params = { .w = w, .t = w->tris[node->i]};
      node->norm_f = t->vol/4.0*integrate( &norm_f_s, &params, -1, 1);
    }
  #else
    tri *t = w->tris[node->i];
    norm_f_s_p params = { .w = w, .t = w->tris[node->i]};
    node->norm_f = t->vol/4.0*integrate( &norm_f_s, &params, -1, 1);
  #endif
  return node->norm_f;
}

/* compute e_n( \node) by computing \|f\|^2_{2,\node} - \|p_n\|^2_{2,\node} */
double error( workspace *w, tree *node, int r) {
  assert( r > 0);
  //find polynomial
  int n = find_n(r);
  #if CACHE_ERROR
    if( node->hp) {
      double e = node->info.hp->e[(n+2)*(n+1)/2-1];
      if( e > -1) {
        return e;
      }
    }
  #endif
  double norm_p_n_sq = construct_p_n( w, node, (n+2)*(n+1)/2);

  find_norm_f( w, node);
  //fprintf( stderr, "error %i,%i (=%i) = %g, norm = %g\n", node->i, r, (n+2)*(n+1)/2, node->norm_f, norm_p_n_sq);
  return node->norm_f - norm_p_n_sq;
}
