#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "tree.h"
#include "partition.h"
#include "edge.h"
#include "tri.h"
#include "workspace.h"

tree *tree_create( int i, tree *parent, tree *left, tree *right) {
  tree *new = malloc( sizeof( tree));
  new->i = i;
  new->parent = parent;

  new->left = left;
  new->right = right;
  return new;
}

int tree_is_leaf( tree *node) {
  return node->left == NULL;
}

void tree_subdivide( workspace *w, tree *node) {
  assert( tree_is_leaf( node));

  //create point midway refinement edge
  tri *t = w->tris[node->i];
  point b = w->points[t->p[1]];
  point c = w->points[t->p[2]];
  point p = {.x = (b.x+c.x)/2, .y = (b.y+c.y)/2};
  int i = workspace_add_point( w, p);

  //create tris
  tri *left = tri_create( w, i, t->p[0], t->p[1]);
  tri *right = tri_create( w, i, t->p[2], t->p[0]);
  int li = workspace_add_tri( w, left);
  int ri = workspace_add_tri( w, right);

  //create trees
  tree *lt = tree_create( li, node, NULL, NULL);
  tree *rt = tree_create( ri, node, NULL, NULL);
  node->left = lt;
  node->right = rt;

  //change leaves
  workspace_remove_leaf( w, node);
  //TODO: we have now removed a leaf, so every element of the edge matrix with index > leaf index of node is now off by one.
  int lli = workspace_add_leaf( w, lt);
  int rli = workspace_add_leaf( w, rt);

  //change edge matrix
  edge_reset( w->edges, t->p[1], t->p[2]);
  edge_set( w->edges, t->p[0], t->p[1], lli);
  edge_set( w->edges, t->p[1], i, lli);
  edge_set( w->edges, i, t->p[0], lli);
  edge_set( w->edges, t->p[2], t->p[0], rli);
  edge_set( w->edges, t->p[0], i, rli);
  edge_set( w->edges, i, t->p[2], rli);
}

void tree_free_shallow( tree *node) {
  free( node);
}

void tree_free_deep( tree *node) {
  if( !tree_is_leaf( node)) {
    tree_free_deep( node->left);
    tree_free_deep( node->right);
  }
  tree_free_shallow( node);
}

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

  for( i = 0; i < w->nleaves; i++) {
    if( w->leaves[i]->i == 12) {
      tree *list[1] = {w->leaves[i]};
      partition_refine( w, list, 1);

      break;
    }
  }

  workspace_print( w);

  /*
  for( i = 0; i < w->nleaves; i++) {
    if( w->leaves[i]->i == 15) {
      tree *list[1] = {w->leaves[i]};
      partition_refine( w, list, 1);

      break;
    }
  }

  workspace_print( w);
  */
}
