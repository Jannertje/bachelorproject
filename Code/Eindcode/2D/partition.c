#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "partition.h"
#include "tree.h"
#include "tri.h"
#include "edge.h"
#include "workspace.h"

/* given tris, create trees and setup initial partition */
void partition_setup( workspace *w, int hp) {
  int i, j;
  for( i = 0; i < w->ntris; i++) {
    //create tree from tri and add root
    tree *node = tree_create( w, i, NULL, NULL, NULL, hp);
    workspace_add_root( w, node);

    //add leaf
    workspace_add_leaf( w, w->roots[i]);

    //set edge matrix
    for( j = 0; j < 3; j++) {
      edge_set( w->edges, w->tris[i]->p[j], w->tris[i]->p[(j+1)%3], w->roots[i]);
    }
  }
}

/* helper function to partition_match() */
void partition_match_node( workspace *w, tree *node) {
  //subdivide
  tree_subdivide( w, node);

  //set the tris correctly
  tri *left = w->tris[node->left->i];
  tri *right = w->tris[node->right->i];
  tri *newleft = tri_create( w, left->p[2], left->p[0], left->p[1]);
  tri *newright = tri_create( w, right->p[1], right->p[2], right->p[0]);
  w->tris[node->left->i] = newleft;
  w->tris[node->right->i] = newright;
  tri_free( left);
  tri_free( right);
}

/* create a matching partition by means of Algorithm A.11 */
void partition_match( workspace *w) {
  int i;
  if( w->is_matched) return;
  for( i = 0; i < w->nroots; i++) {
    tree_subdivide( w, w->roots[i]);
    partition_match_node( w, w->roots[i]->left);
    partition_match_node( w, w->roots[i]->right);
  }
  //substitute all roots by the newly made partition
  w->nroots = 0;
  for( i = 0; i < w->nleaves; i++) {
    workspace_add_root( w, w->leaves[i]);
  }

  w->is_matched = 1;
  w->is_conform = 1;
}

/* Create a conforming refinement where node is refined. Algorithm A.5 */
void partition_refine_node( workspace *w, tree *node) {
  assert( w->is_conform);
  assert( tree_is_leaf( node));
  tri *t = w->tris[node->i];
  tree *neighbour = edge_get( w->edges, t->p[2], t->p[1]);
  if( neighbour == NULL) {
    //no neighbouring triangle; we are on the edge
    tree_subdivide( w, node);
  } else {
    tri *neighbour_tri = w->tris[neighbour->i];
    if( neighbour_tri->p[2] == t->p[1] && neighbour_tri->p[1] == t->p[2]) {
      //same refinement edge
      tree_subdivide( w, node);
      tree_subdivide( w, neighbour);
    } else {
      partition_refine_node( w, neighbour);
      tree *left = neighbour->left;
      tree *right = neighbour->right;
      tree *children[2] = {left, right};
      int i;
      for( i = 0; i < 2; i++) {
        if( w->tris[children[i]->i]->p[1] == t->p[2] && 
            w->tris[children[i]->i]->p[2] == t->p[1]) {
          //we found the right child
          tree_subdivide( w, children[i]);
          tree_subdivide( w, node);
        }
      }
    }
  }

  assert( !tree_is_leaf( node));

  w->is_conform = 1;
}

/* Create a conforming refinement where nodes[n] are refined. Algorithm A.13 */
void partition_refine( workspace *w, tree **nodes, int n) {
  assert( w->is_conform);
  int i;
  for( i = 0; i < n; i++) {
    if( tree_is_leaf( nodes[i])) {
      partition_refine_node( w, nodes[i]);
    }
  }
}

/* Algorithm A.7 */
void partition_make_conform( workspace *w) {
  if( w->is_conform) return;

  int i, j, k = 0, len = 0;
  tree **to_refine = NULL;
  //create list to refine
  for( i = 0; i < w->nleaves; i++) {
    tree *leaf = w->leaves[i];
    if( tree_has_hanging_vertex( w, leaf)) {
      if( len == 0) {
        to_refine = malloc( 2 * sizeof( tree *));
        len = 2;
      }

      if( len <= k) {
        to_refine = realloc( to_refine, 2 * len * sizeof( tree *));
        len *= 2;
      }

      to_refine[k++] = leaf;
    }
  }

  //traverse
  while( k > 0) {
    tree *node = to_refine[0];
    for( j = 0; j < k-1; j++) {
      to_refine[j] = to_refine[j+1];
    }
    k--;
    tree_subdivide( w, node);

    //carry over degree
    if( node->hp) {
      node->left->info.hp->r = node->info.hp->r;
      node->right->info.hp->r = node->info.hp->r;
      hptree_get_e( w, node->left, node->info.hp->r);
      hptree_get_e( w, node->right, node->info.hp->r);
    } else {
      node->left->info.h->r = node->info.h->r;
      node->right->info.h->r = node->info.h->r;
    }

    tri *nt = w->tris[node->i];
    if( tree_has_hanging_vertex( w, node->left)) {
      if( len <= k) {
        to_refine = realloc( to_refine, 2 * len * sizeof( tree *));
        len *= 2;
      }

      to_refine[k++] = node->left;
    }
    if( tree_has_hanging_vertex( w, node->right)) {
      if( len <= k) {
        to_refine = realloc( to_refine, 2 * len * sizeof( tree *));
        len *= 2;
      }

      to_refine[k++] = node->right;
    }

    //last step
    for( i = 0; i < w->nleaves; i++) {
      if( w->leaves[i] == node) continue;

      //leaf in list already?
      int found = 0;
      for( j = 0; j < k; j++) {
        if( to_refine[j] == w->leaves[i]) found = 1;
      }
      if( found) continue;

      tri *t = w->tris[w->leaves[i]->i];
      if( edge_get( w->edges, t->p[1], t->p[0]) == node ||
          edge_get( w->edges, t->p[2], t->p[1]) == node ||
          edge_get( w->edges, t->p[0], t->p[2]) == node ||
          edge_get( w->edges, nt->p[1], nt->p[0]) == w->leaves[i] ||
          edge_get( w->edges, nt->p[2], nt->p[1]) == w->leaves[i] ||
          edge_get( w->edges, nt->p[0], nt->p[2]) == w->leaves[i]) {
        //we are neighbours
        if( tree_has_hanging_vertex( w, w->leaves[i])) {
          if( len <= k) {
            to_refine = realloc( to_refine, 2 * len * sizeof( tree *));
            len *= 2;
          }

          to_refine[k++] = w->leaves[i];
        }
      }
    }
  }
  w->is_conform = 1;
}

/* Find inner nodes of the total partition. */
void partition_inner_nodes( workspace *w, tree ***inners, int *ninners) {
  tree **queue = malloc( w->nleaves * sizeof( tree *));
  *inners = malloc( w->nleaves * sizeof( tree *));
  *ninners = 0;
  int nqueue = w->nroots;

  int i;
  for( i = 0; i < w->nroots; i++) {
    queue[i] = w->roots[i];
  }

  int j;
  while( nqueue > 0) {
    tree *node = queue[0];
    for( j = 0; j < nqueue - 1; j++) {
      queue[j] = queue[j+1];
    }
    nqueue--;

    if( !tree_is_leaf( node)) {
      for( j = *ninners; j >= 1; j--) {
        (*inners)[j] = (*inners)[j-1];
      }
      (*inners)[0] = node;
      (*ninners)++;
      queue[nqueue++] = node->left;
      queue[nqueue++] = node->right;
    }
  }

  free( queue);
}
