#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "tree.h"
#include "helper.h"
#include "error.h"
#include "partition.h"
#include "edge.h"
#include "tri.h"
#include "treegen.h"
#include "workspace.h"

hp *hp_create() {
  hp *new = malloc( sizeof( hp));
  new->r = -1;
  new->te = -1;
  new->q = -1;

  new->ehp = malloc( 2*sizeof( double));
  new->ehp[0] = -1;
  new->ehp[1] = -1;
  new->lenehp = 2;

  new->tehp = malloc( 2*sizeof( double));
  new->tehp[0] = -1;
  new->tehp[1] = -1;
  new->lentehp = 2;

  new->e = malloc( 2*sizeof( double));
  new->e[0] = -1;
  new->e[1] = -1;
  new->lene = 2;

  return new;
}

tree *tree_create( workspace *w, int i, tree *parent, tree *left, tree *right, int hp) {
  tree *new = malloc( sizeof( tree));
  new->i = i;
  new->parent = parent;

  new->left = left;
  new->right = right;
  new->hp = hp;
  new->t = NULL;
  new->gammas = malloc( 2*sizeof( double));
  new->gammas[0] = NAN;
  new->gammas[1] = NAN;
  /*
  new->gammas[2] = NAN;
  new->gammas[3] = NAN;
  new->gammas[4] = NAN;
  new->gammas[5] = NAN;
  new->gammas[6] = NAN;
  new->gammas[7] = NAN;
  new->gammas[8] = NAN;
  new->gammas[9] = NAN;
  */
  new->lengammas = 2;
  new->hgammas = 0;

  if( hp) {
    new->info.hp = hp_create();
  } else {
    new->info.h = malloc( sizeof( h));
    new->info.h->r = -1;
    new->info.h->e = -1;
    new->info.h->te = -1;
  }

  return new;
}

double hptree_get_ehp( workspace *w, tree *node, int r) {
  assert( node->hp);
  if( node->info.hp->lenehp < r || node->info.hp->ehp[r-1] == -1) {
    printf("watwat ehp te klein\n");
    exit(1);
  }
  return node->info.hp->ehp[r-1];
}

void hptree_set_ehp( workspace *w, tree *node, int r, double val) {
  assert( node->hp);

  if( node->info.hp->lenehp < r) {
    int i;
    int pow2 = pow2roundup( r);
    node->info.hp->ehp = realloc( node->info.hp->ehp, 2*pow2*sizeof( double));
    for( i = node->info.hp->lenehp; i < 2*pow2; i++) {
      node->info.hp->ehp[i] = -1;
    }
    node->info.hp->lenehp = 2*pow2;
  }
  node->info.hp->ehp[r-1] = val;
}

double hptree_get_tehp( workspace *w, tree *node, int r) {
  assert( node->hp);
  if( node->info.hp->lentehp < r || node->info.hp->tehp[r-1] == -1) {
    printf("watwat tehp te klein\n");
    exit(1);
  }
  return node->info.hp->tehp[r-1];
}

void hptree_set_tehp( workspace *w, tree *node, int r, double val) {
  assert( node->hp);

  if( node->info.hp->lentehp < r) {
    int i;
    int pow2 = pow2roundup( r);
    node->info.hp->tehp = realloc( node->info.hp->tehp, 2*pow2*sizeof( double));
    for( i = node->info.hp->lentehp; i < 2*pow2; i++) {
      node->info.hp->tehp[i] = -1;
    }
    node->info.hp->lentehp = 2*pow2;
  }
  node->info.hp->tehp[r-1] = val;
}

double hptree_get_e( workspace *w, tree *node, int r) {
  assert( node->hp);

  if( node->info.hp->lene < r) {
    int i;
    int pow2 = pow2roundup( r);
    node->info.hp->e = realloc( node->info.hp->e, 2*pow2*sizeof( double));
    for( i = node->info.hp->lene; i < 2*pow2; i++) {
      node->info.hp->e[i] = -1;
    }
    node->info.hp->lene = 2*pow2;
  }
  if( node->info.hp->e[r-1] > -1) {
    return node->info.hp->e[r-1];
  } else {
    return (node->info.hp->e[r-1] = error( w, node, r));
  }
}

double htree_get_e( workspace *w, tree *node) {
  assert( !node->hp);
  if( node->info.h->e > -1) {
    return node->info.h->e;
  } else {
    return (node->info.h->e = error( w, node, node->info.h->r));
  }
}

void tree_set_gamma( workspace *w, tree *node, int i, double gamma) {
  assert( node->lengammas > i);

  node->gammas[i] = gamma;
  if( i > node->hgammas) {
    node->hgammas = i;
  }
}

double tree_get_gamma( workspace *w, tree *node, int i) {
  if( node->lengammas <= i) {
    int j;
    int pow2 = pow2roundup( i);
    node->gammas = realloc( node->gammas, 2*pow2*sizeof( double));
    for( j = node->hgammas+1; j < 2*pow2; j++) {
      node->gammas[j] = NAN;
    }
    node->lengammas = 2*pow2;
  }

  return node->gammas[i];
}

int tree_is_leaf( tree *node) {
  return node->left == NULL;
}

void tree_trim( workspace *w, tree *node) {
  assert( !tree_is_leaf( node));
  if( !tree_is_leaf( node->left)) {
    tree_trim( w, node->left);
  } else {
    workspace_remove_leaf( w, node->left);
  }
  node->left = NULL;
  if( !tree_is_leaf( node->right)) {
    tree_trim( w, node->right);
  } else {
    workspace_remove_leaf( w, node->right);
  }
  node->right = NULL;
  workspace_add_leaf( w, node);
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
  tree *lt = tree_create( w, li, node, NULL, NULL, node->hp);
  tree *rt = tree_create( w, ri, node, NULL, NULL, node->hp);
  node->left = lt;
  node->right = rt;

  //change leaves
  workspace_remove_leaf( w, node);
  workspace_add_leaf( w, lt);
  workspace_add_leaf( w, rt);

  //change edge matrix
  edge_set( w->edges, t->p[0], t->p[1], lt);
  edge_set( w->edges, i, t->p[0], lt);
  edge_set( w->edges, t->p[2], t->p[0], rt);
  edge_set( w->edges, t->p[0], i, rt);

  //set error
  if( !node->hp && node->info.h->r > 0) {
    lt->info.h->r = node->info.h->r;
    lt->info.h->te = inverror( htree_get_e( w, lt), node->info.h->te);
    rt->info.h->r = node->info.h->r;
    rt->info.h->te = inverror( htree_get_e( w, rt), node->info.h->te);
  }

  w->is_conform = 0;
}

void htree_sort_insertion( tree **nodes, int n) {
  int i, k;
  tree *temp;
  for( i = 1; i < n; i++) {
    for( k = i; k > 0 && nodes[k]->info.h->te > nodes[k-1]->info.h->te; k--) {
      temp = nodes[k];
      nodes[k] = nodes[k-1];
      nodes[k-1] = temp;
    }
  }
}

int tree_is_on_edge( workspace *w, tree *node) {
  int j;
  for( j = 0; j < 3; j++) {
    tree *neighbour = edge_get( 
      w->edges, 
      w->tris[node->i]->p[(j+1)%3], 
      w->tris[node->i]->p[j]
    );
    if( neighbour == NULL) {
      //we are on the edge
      return 1;
    }
  }
  return 0;
}

int tree_has_hanging_vertex( workspace *w, tree *node) {
  int j;
  for( j = 0; j < 3; j++) {
    tree *neighbour = edge_get( 
      w->edges, 
      w->tris[node->i]->p[(j+1)%3], 
      w->tris[node->i]->p[j]
    );
    if( neighbour != NULL && 
        !tree_is_leaf( neighbour) && 
        (!tree_is_leaf( neighbour->left) || 
         (w->tris[neighbour->i]->p[0] != w->tris[node->i]->p[j] &&
          w->tris[neighbour->i]->p[0] != w->tris[node->i]->p[(j+1)%3]))) {
      //see stevensons MakeConform
      return 1;
    }
  }
  return 0;
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

  for( i = 0; i < w->nleaves; i++) {
    if( w->leaves[i]->i == 12) {
      tree *list[1] = {w->leaves[i]};
      partition_refine( w, list, 1);

      break;
    }
  }


  tree *list[2] = {NULL};
  for( i = 0; i < w->nleaves; i++) {
    if( w->leaves[i]->i == 15) {
      list[0] = w->leaves[i];
    } else if( w->leaves[i]->i == 17) {
      list[1] = w->leaves[i];
    }
  }
  partition_refine( w, list, 2);
  workspace_print_plot( w);
}
*/
