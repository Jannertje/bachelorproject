#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "partition.h"
#include "tree.h"
#include "tri.h"
#include "edge.h"
#include "workspace.h"

void partition_match( workspace *w);

void partition_setup( workspace *w) {
  int i, j;
  for( i = 0; i < w->ntris; i++) {
    //create tree from tri and add root
    workspace_add_root( w, tree_create( i, NULL, NULL, NULL));
    //add leaf
    workspace_add_leaf( w, w->roots[i]);

    //set edge matrix
    for( j = 0; j < 3; j++) {
      edge_set( w->edges, w->tris[i]->p[j], w->tris[i]->p[(j+1)%3], i);
    }
  }

  //make a matching partition
  partition_match( w);
}

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

void partition_match( workspace *w) {
  int i;
  for( i = 0; i < w->nroots; i++) {
    tree_subdivide( w, w->roots[i]);
    partition_match_node( w, w->roots[i]->left);
    partition_match_node( w, w->roots[i]->right);
  }
}

void partition_refine_node( workspace *w, tree *node) {
  workspace_print( w);
  printf("refining node with tri number %i\n", node->i);
  tri *t = w->tris[node->i];
  int neighbour_leaf_index = edge_get_leaf( w->edges, t->p[2], t->p[1]);
  printf("neighbour tri index %i\n", w->leaves[neighbour_leaf_index]->i);
  if( neighbour_leaf_index == -1) {
    //no neighbouring triangle; we are on the edge
    printf("we are on the edge\n");
    tree_subdivide( w, node);
    return;
  } else {
    tree *neighbour = w->leaves[neighbour_leaf_index];
    tri *neighbour_tri = w->tris[neighbour->i];
    if( neighbour_tri->p[2] == t->p[1] && neighbour_tri->p[1] == t->p[2]) {
      //same refinement edge
      printf("%i and %i have same refinement edge\n", neighbour->i, node->i);
      tree_subdivide( w, node);
      tree_subdivide( w, neighbour);
      return;
    } else {
      printf("entering recursion..\n");
      partition_refine_node( w, neighbour);
      printf("exit recursion.\n");
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
          return;
        }
      }
    }
  }
  printf("watwat\n");
  exit(1);
}

void partition_refine( workspace *w, tree **nodes, int n) {
  int i;
  for( i = 0; i < n; i++) {
    if( tree_is_leaf( nodes[i])) {
      partition_refine_node( w, nodes[i]);
    }
  }
}
