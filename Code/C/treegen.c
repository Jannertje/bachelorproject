#include <stdio.h>
#include <gsl/gsl_math.h>
#include <stdlib.h>
#include <assert.h>
#include "types.h"
#include "helper.h"
#include "treegen.h"
#include "tree.h"
#include "partition.h"
#include "workspace.h"
#include "error.h"
#include "math.h"

double inverror( double a, double b) {
  if( a == 0.0 || b == 0.0) return 0.0;
  else return 1/(1/a + 1/b);
}

void treegen_hp( workspace *w, int Nmax) {
  int i, j, n = 0, len = 2;
  tree **to_refine = malloc( 2 * sizeof( tree *));
  for( i = 0; i < w->nroots; i++) {
    int N = 1;
    //step 1
    //printf("step 1\n");
    tree *D = w->roots[i];
    assert( D->hp);
    D->info.hp->r = 1;
    D->info.hp->te = hptree_get_e( w, D, 1);
    hptree_set_ehp( w, D, 1, D->info.hp->te);
    hptree_set_tehp( w, D, 1, D->info.hp->te);
    D->info.hp->q = D->info.hp->te;
    D->t = D;

    outerloop: {
      //step 2
      //printf("step 2\n");
      tree *node_N = D->t;
      //printf("  node_N := %i\n", node_N->i);
      tree_subdivide( w, node_N);
      //printf("  subdivision yielded %i,%i\n", node_N->left->i, node_N->right->i);
      if( n + 1 >= len) {
        int pow2 = pow2roundup( len);
        to_refine = realloc( to_refine, 2*pow2*sizeof( tree *));
        len *= 2;
      }
      to_refine[n++] = node_N->left;
      to_refine[n++] = node_N->right;

      //step 3
      //printf("step 3\n");
      tree *children[2] = {node_N->left, node_N->right};
      for( j = 0; j < 2; j++) {
        tree *node = children[j];
        node->info.hp->r = 1;
        node->info.hp->te = inverror( hptree_get_e( w, node, 1), node->info.hp->te);
        hptree_set_ehp( w, node, 1, hptree_get_e( w, node, 1));
        hptree_set_tehp( w, node, 1, node->info.hp->te);
        node->info.hp->q = node->info.hp->te;
        node->t = node;
      }

      //step 4
      //printf("step 4\n");
      tree *node = node_N;

      //step 5
      //printf("step 5\n");
      N += 1;
      if( N >= Nmax) {
        continue;
      }

      //step 6
      //printf("step 6\n");
      innerloop: {
        //printf("  node = %i\n", node->i);
        //step 6a
        node->info.hp->r += 1;
        //printf(" upping r to %i\n", node->info.hp->r);
        
        //step 6b
        tree *node_1 = node->left;
        tree *node_2 = node->right;

        //step 6c
        hptree_set_ehp( 
          w, node, node->info.hp->r, 
          fmin( 
            hptree_get_ehp( w, node_1, node_1->info.hp->r) + 
            hptree_get_ehp( w, node_2, node_2->info.hp->r),
            hptree_get_e( w, node, node->info.hp->r)
          )
        );

        //step 6d
        hptree_set_tehp( 
          w, node, node->info.hp->r, 
          inverror( 
            hptree_get_ehp( w, node, node->info.hp->r),
            hptree_get_ehp( w, node, node->info.hp->r-1)
          )
        );

        //step 6e
        tree *X = node_2;
        if( node_1->info.hp->q > node_2->info.hp->q) X = node_2;
        node->info.hp->q = fmin( X->info.hp->q, hptree_get_tehp( w, node, node->info.hp->r));
        node->t = X->t;

        //step 6f
        if( node == D) {
          goto outerloop;
        } else {
          node = node->parent;
          goto innerloop;
        }
      }
    }
  }

  //partition_refine( w, to_refine, n);
  free( to_refine);

  //create T_N from \T_N.
  tree **inners = NULL;
  int ninners;
  partition_inner_nodes( w, &inners, &ninners);
  for( i = 0; i < ninners; i++) {
    double ehp = hptree_get_ehp( w, inners[i], inners[i]->info.hp->r);
    double er = hptree_get_e( w, inners[i], inners[i]->info.hp->r);
    if( !gsl_fcmp( ehp, er, 1E-5)) {
      tree_trim( w, inners[i]);
    }
  }
}

void treegen_h( workspace *w, int Nmax, int r) {
  int i, N;
  for( i = 0; i < w->nroots; i++) {
    tree *D = w->roots[i];
    D->info.h->r = r;
    D->info.h->te = htree_get_e( w, D);
  }
  htree_sort_insertion( w->roots, w->nroots);

  for( N = 0; N < Nmax; N++) {
    int i, j = 1;
    double te = w->leaves[0]->info.h->te;
    for( i = 0; i < w->nleaves; i++) {
      if( w->leaves[i]->info.h->te != te) {
        j = i+1;
        break;
      }
    }
    for( i = 0; i < j; i++) {
      tree *parent = w->leaves[i];
      tree_subdivide( w, parent);
    }
    htree_sort_insertion( w->leaves, w->nleaves);
  }
}
