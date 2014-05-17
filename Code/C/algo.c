#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "algo.h"
#include "tree.h"
#include "error.h"

typedef enum algo_type { GREEDY, BINEV2004, BINEV2007 } algo_type;
typedef error_info ( *algo_fun)( algo_info *i, tree *node, int r);
tree *runner( algo_type algo, algo_info *i, int r, int n);
error_info greedy_error( algo_info *i, tree *node, int r);
error_info binev2004_error( algo_info *i, tree *node, int r);
error_info binev2007_error( algo_info *i, tree *node, int r);

error_info greedy_error( algo_info *i, tree *node, int r) {
  if( node->e.error > -1.0) {
    return node->e;
  } else {
    double error = i->e( i->f, node->b, r);
    error_info e = {error, error};
    node->e = e;
    return e;
  }
}

tree *algo_greedy( algo_info *i, int r, int n) {
  return runner( GREEDY, i, r, n);
}

/*
 * Wrapper around the error function from error.h.
 */
double real_error( algo_info *i, tree *node, int r) {
  if( node->e.real_error > -1.0) {
    return node->e.real_error;
  } else {
    double real_error = i->e( i->f, node->b, r);
    node->e.real_error = real_error;
    return real_error;
  }
}

double binev2004_q( algo_info *i, tree *node, int r) {
  assert( !tree_is_leaf( node));

  double teller = real_error( i, node->forest[0], r) 
                + real_error( i, node->forest[1], r);
  error_info error = binev2004_error( i, node, r);
  return teller/(error.error + error.real_error)*error.error;
}

error_info binev2004_error( algo_info *i, tree *node, int r) {
  if( node->e.error > -1.0) {
    return node->e;
  } else {
    double realerror = real_error( i, node, r);
    double error = realerror;
    if( !tree_is_root( node)) {
      error = binev2004_q( i, node->parent, r);
    }
    error_info e = {.error = error, 
                    .real_error = realerror};
    node->e = e;
    return e;
  }
}

tree *algo_binev2004( algo_info *i, int r, int n) {
  return runner( BINEV2004, i, r, n);
}

error_info binev2007_error( algo_info *i, tree *node, int r) {
  if( node->e.error > -1.0) {
    return node->e;
  } else {
    double realerror = real_error( i, node, r);
    double error = realerror;
    if( !tree_is_root( node)) {
      error_info parenterror = binev2007_error( i, node->parent, r);
      error = 1/( 1/realerror + 1/parenterror.error);
    }
    error_info e = {.error = error, .real_error = realerror};
    node->e = e;
    return e;
  }
}

tree *algo_binev2007( algo_info *i, int r, int n) {
  return runner( BINEV2007, i, r, n);
}

int iterator( algo_fun errorfun, algo_info *i, int r, tree **t) {
  tree_list *leaves = tree_leaves( *t);
  tree_list *cur = leaves;

  while( cur != NULL) {
    cur->node->e = errorfun( i, cur->node, r);
    cur = cur->next;
  }

  tree_list *bests = NULL;
  i->s( leaves, &bests);

  if( bests == NULL) {
    return 1;
  }

  cur = bests;
  while( cur != NULL) {
    tree_subdivide( cur->node);
    cur = cur->next;
  }

  tree_list_free( leaves);
  tree_list_free( bests);
  return 0;
}

tree *runner( algo_type algo, algo_info *i, int r, int n) {
  algo_fun errorfun;
  switch( algo) {
    case GREEDY:
      errorfun = greedy_error;
      printf("Running Greedy algorithm, degrees of freedom %i; iterations %i.\n", r, n);
      break;
    case BINEV2004:
      errorfun = binev2004_error;
      printf("Running Binev2004 algorithm, degrees of freedom %i; iterations %i.\n", r, n);
      break;
    case BINEV2007:
      errorfun = binev2007_error;
      printf("Running Binev2007 algorithm, degrees of freedom %i; iterations %i.\n", r, n);
      break;
  }
  
  error_info e = {-1.0, -1.0};
  tree *t = tree_create( i->b, NULL, NULL, NULL, e);

  if( n < 0) {
    int j = 0;
    while( 1) {
      printf("Iteration %i:\n", j);
      if( iterator( errorfun, i, r, &t)) {
        break;
      }
      j++;
    }
  } else {
    for( int j = 0; j < n; j++) {
      printf("Iteration %i:\n", j);
      if( iterator( errorfun, i, r, &t)) {
        break;
      }
    }
  }

  double leaves_sum_real_error = 0.0;
  tree_list *leaves = tree_leaves( t);
  tree_list *cur = leaves;
  while( cur != NULL) {
    if( cur->node->e.real_error == -1.0) {
      cur->node->e.real_error = i->e( i->f, cur->node->b, r);
    }
    leaves_sum_real_error += cur->node->e.real_error;
    cur = cur->next;
  }
  tree_list_free( leaves);
  tree_print( t);

  printf("Total sum of leaf errors: %g\n", leaves_sum_real_error);

  return t;
}

