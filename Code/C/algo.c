#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "algo.h"
#include "tree.h"
#include "error.h"

typedef enum algo_type { GREEDY, BINEV2004, BINEV2007 } algo_type;
typedef error_info ( *algo_fun)( algo_info *i, tree *node, int r);
tree *runner( algo_type algo, algo_info *i, int r, int n);
error_info greedy_error( algo_info *i, tree *node, int r);
error_info binev2004_error( algo_info *i, tree *node, int r);
error_info binev2007_error( algo_info *i, tree *node, int r);

/* ------------------- COMMON ----------------- */

/*
 * Wrapper around the error function from error.h.
 */
double real_error( algo_info *i, tree *node, int r) {
  double e;
  if( (e = htree_error_info( node).real_error) > -1.0) {
    return e;
  } else {
    double real_error = i->e( i->f, i->b, node->l, r);
    node->t.h->e.real_error = real_error;
    return real_error;
  }
}

/*
 * Do one iteration of the $h$-type algorithms. They are so alike we need just 1
 * function.
 */
int iterator( algo_fun errorfun, algo_info *i, int r, tree **t) {
  tree_list *leaves = tree_leaves( *t);
  tree_list *cur = leaves;

  while( cur != NULL) {
    cur->node->t.h->e = errorfun( i, cur->node, r);
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
      printf("Running Greedy algorithm, degrees of "
             "freedom %i; iterations %i.\n", r, n);
      break;
    case BINEV2004:
      errorfun = binev2004_error;
      printf("Running Binev2004 algorithm, degrees of "
             "freedom %i; iterations %i.\n", r, n);
      break;
    case BINEV2007:
      errorfun = binev2007_error;
      printf("Running Binev2007 algorithm, degrees of "
             "freedom %i; iterations %i.\n", r, n);
      break;
  }
  
  error_info e = {-1.0, -1.0};
  location init = {0, 0};
  tree *t = htree_create( init, NULL, NULL, NULL, e);

  if( n < 0) {
    int j = 0;
    while( 1) {
      if( iterator( errorfun, i, r, &t)) {
        break;
      }
      j++;
    }
  } else {
    int j;
    for( j = 0; j < n; j++) {
      if( iterator( errorfun, i, r, &t)) {
        break;
      }
    }
  }

  double leaves_sum_real_error = 0.0;
  tree_list *leaves = tree_leaves( t);
  tree_list *cur = leaves;
  while( cur != NULL) {
    if( htree_error_info( cur->node).real_error == -1.0) {
      cur->node->t.h->e.real_error = i->e( i->f, i->b, cur->node->l, r);
    }
    leaves_sum_real_error += cur->node->t.h->e.real_error;
    cur = cur->next;
  }
  tree_list_free( leaves);
  tree_print( t);

  printf("Total sum of leaf errors: %g\n", leaves_sum_real_error);

  return t;
}


/* ------------------- GREEDY ---------------- */

/* 
 * Construct the greedy error_info.error == error_info.real_error
 */
error_info greedy_error( algo_info *i, tree *node, int r) {
  error_info e;
  if( (e = htree_error_info( node)).error > -1.0) {
    return e;
  } else {
    double error = i->e( i->f, i->b, node->l, r);
    e.error = error;
    e.real_error = error;
    node->t.h->e = e;
    return e;
  }
}

/*
 * Wrapper around runner() to allow controlled call by externals
 */
tree *algo_greedy( algo_info *i, int r, int n) {
  return runner( GREEDY, i, r, n);
}

/* ------------------- BINEV2004 ---------------- */

/*
 * See binev2004 algorithm in thesis.
 */
double binev2004_q( algo_info *i, tree *node, int r) {
  assert( !tree_is_leaf( node));

  double teller = real_error( i, node->forest[0], r) 
                + real_error( i, node->forest[1], r);
  error_info error = binev2004_error( i, node, r);
  return teller/(error.error + error.real_error)*error.error;
}

/*
 * Construct the binev2004 modified error. See thesis.
 */
error_info binev2004_error( algo_info *i, tree *node, int r) {
  error_info e;
  if( (e = htree_error_info( node)).error > -1.0) {
    return e;
  } else {
    double realerror = real_error( i, node, r);
    double error = realerror;
    if( !tree_is_root( node)) {
      error = binev2004_q( i, node->parent, r);
    }
    e.error = error;
    e.real_error = realerror;
    node->t.h->e = e;
    return e;
  }
}

/*
 * Wrapper around runner().
 */
tree *algo_binev2004( algo_info *i, int r, int n) {
  return runner( BINEV2004, i, r, n);
}

/* --------------- BINEV2007 -------------------- */

/* 
 * Construct Binev2007 modified error. See thesis.
 */
error_info binev2007_error( algo_info *i, tree *node, int r) {
  error_info e;
  if( (e = htree_error_info( node)).error > -1.0) {
    return e;
  } else {
    double realerror = real_error( i, node, r);
    double error = realerror;
    if( !tree_is_root( node)) {
      error_info parenterror = binev2007_error( i, node->parent, r);
      error = 1/( 1/realerror + 1/parenterror.error);
    }
    e.error = error;
    e.real_error = realerror;
    node->t.h->e = e;
    return e;
  }
}

tree *algo_binev2007( algo_info *i, int r, int n) {
  return runner( BINEV2007, i, r, n);
}

/* --------------------- BINEV2013 ----------------- */

/*
 * Find Binev2013 initial modified error (step 2 of algorithm)
 */
error_info binev2013_error( algo_info *i, tree *node) {
  error_info e;
  if( (e = htree_error_info( node)).error > -1.0) {
    return e;
  } else {
    double realerror = real_error( i, node, 1);
    double error = realerror;
    if( !tree_is_root( node)) {
      error_info parenterror = binev2007_error( i, node->parent, 1);
      error = 1/( 1/realerror + 1/parenterror.error);
    }
    e.error = error;
    e.real_error = realerror;
    node->t.h->e = e;
    return e;
  }
}

/*
 * Find total error of a node. See Thesis for lemmata to prove this works.
 */
double binev2013_total_error( algo_info *i, tree *node, tree *T, tree *TN) {
  tree *TN_node = tree_find_node( node, TN);
  assert( TN_node != NULL);
  //find elements in X = TN_node \cap T as subset of TN
  tree_list *X = tree_list_create( TN_node); //we know that TN_node \in X
  tree_list *LX = NULL;
  tree_list *cur = X;

  //we want to traverse all nodes in X
  tree_list *queue = tree_list_create( TN_node);
  while( queue != NULL) {
    tree *queue_node = tree_list_popleft( &queue); //queue_node is in X
    if( !tree_is_leaf( queue_node) && 
        tree_find_node( queue_node->forest[0], T) != NULL) { //queue_node is in I(X)
      //queue_node->forest[0] is in T
      if( queue == NULL) {
        queue = tree_list_create( queue_node->forest[0]);
      } else {
        //TODO: keep a curqueue var for O(1) insert
        tree_list_append( queue, queue_node->forest[0]);
      }
      //if forest[0] is in T, then so is forest[1]
      tree_list_append( queue, queue_node->forest[1]);
      cur = tree_list_insert_after( cur, queue_node->forest[0]);
      cur = tree_list_insert_after( cur, queue_node->forest[1]); 
    } else {
      if( LX == NULL) {
        LX = tree_list_create( queue_node);
      } else {
        LX = tree_list_insert_beginning( LX, queue_node); //O(1) insert
      }
    }
  }

  double sum = 0;
  cur = LX;
  while( cur != NULL) {
    tree *dummy;
    int r_node = htree_r_node( cur->node, TN, &dummy);
    sum += (i->e)( i->f, i->b, cur->node->l, r_node);
    cur = cur->next;
  }

  tree_list_free( LX);
  tree_list_free( X);
  return sum;
}

/*
 * Step 4 of Binev2013.
 */
tree *binev2013_generator( tree *TN, tree *T, tree *parent) {
  tree *Thp;
  if( tree_is_leaf( T)) {
    tree *dummy = NULL;
    Thp = hptree_create( T->l, NULL, NULL, parent, T->t.h->e, -1);
    Thp->t.hp->r = htree_r_node( T, TN, &dummy);
  } else {
    //no leaf, we dont care about error or r.
    Thp = hptree_create( T->l, NULL, NULL, parent, T->t.h->e, -1);
    Thp->forest[0] = binev2013_generator( TN, T->forest[0], Thp);
    Thp->forest[1] = binev2013_generator( TN, T->forest[1], Thp);
  }
  return Thp;
}

/*
 * Wrapper around recursive function binev2013_generator() to supply no 
 * non-obvious parameters.
 */
tree *binev2013_generate_Thp( tree *TN, tree *T) {
  return binev2013_generator( TN, T, NULL);
}

/* 
 * Full Binev2013 algorithm -- one iteration only.
 * TN is already allocated and Thp will be allocated if generate_Thp is set.
 */
int binev2013_iterator( algo_info *i, tree **TN, int generate_Thp, tree **Thp) {
  assert( !(*TN)->is_hp);

  //step 1
  tree *T = htree_copy( *TN);

  //step 2
  tree_list *leaves = tree_leaves( T);
  tree_list *cur = leaves;

  while( cur != NULL) {
    tree *curnode = tree_find_node( cur->node, *TN);
    curnode->t.h->e = binev2013_error( i, curnode);
    cur = cur->next;
  }

  tree_list_free( leaves);

  //step 3
  tree_list *inners = tree_inner_nodes( T);
  cur = inners;

  while( cur != NULL) {
    tree *TN_node_in_TN = NULL;
    int r_node = htree_r_node( cur->node, *TN, &TN_node_in_TN); 
    //side effect: TN_node_in_TN is now set

    double realerror = (i->e)( i->f, i->b, cur->node->l, r_node);
    double totalerror = binev2013_total_error( i, cur->node, T, *TN);
    if( realerror < totalerror) {
      //step 3a
      tree_trim_subtree( cur->node);
      cur->node->t.h->e.real_error = realerror;

      //step 3b
      tree_list *L_TN_node = tree_leaves( TN_node_in_TN);
      tree_list *curleaf = L_TN_node;
      while( curleaf != NULL) {
        curleaf->node->t.h->e.error *= realerror/totalerror;
        curleaf = curleaf->next;
      }
      tree_list_free( L_TN_node);
    }
    cur = cur->next;
  }

  tree_list_free( inners);

  //step 4
  if( generate_Thp) {
    *Thp = binev2013_generate_Thp( *TN, T);
  }
  tree_free_subtree( T);

  //step 5
  leaves = tree_leaves( *TN);
  cur = leaves;

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

/*
 * Function as called by externals.
 */
tree *algo_binev2013( algo_info *i, int n) {
  printf("Running Binev2013 algorithm; iterations: %i.\n", n);
  error_info e = {-1.0, -1.0};
  location init = {0, 0};
  tree *TN = htree_create( init, NULL, NULL, NULL, e);
  tree *Thp = NULL;

  int generate_Thp = 1; //set this to 1 for the time being.
  if( n < 0) {
    int j = 0;
    while( 1) {
      printf("Iteration %i\n", j);
      if( Thp != NULL) {
        tree_free_subtree( Thp);
      }
      if( binev2013_iterator( i, &TN, generate_Thp, &Thp)) {
        break;
      }
    }
  } else if( n == 0) {
    Thp = binev2013_generate_Thp( TN, TN);
  } else {
    int j;
    for( j = 0; j < n; j++) {
      printf("Iteration %i\n", j);
      if( Thp != NULL) {
        tree_free_subtree( Thp);
      }
      if( binev2013_iterator( i, &TN, generate_Thp, &Thp)) {
        break;
      }
    }
  }

  double leaves_sum_real_error = 0.0;
  tree_list *leaves = tree_leaves( Thp);
  tree_list *cur = leaves;
  while( cur != NULL) {
    if( hptree_error_info( cur->node).real_error < 0.0) {
      cur->node->t.hp->e.real_error = i->e(i->f, 
                                           i->b, 
                                           cur->node->l, 
                                           hptree_r( cur->node));
    }
    leaves_sum_real_error += cur->node->t.hp->e.real_error;
    cur = cur->next;
  }
  tree_list_free( leaves);
  //tree_print( Thp);
  tree_free_subtree( TN);

  printf("Total sum of leaf errors: %g\n", leaves_sum_real_error);

  return Thp;
}
