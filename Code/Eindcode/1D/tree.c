#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "tree.h"
#include "helper.h"
#include "error.h"
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

tree *tree_create( double a, double b, tree *parent, tree *left, tree *right, int hp) {
  assert( b > a);
  tree *new = malloc( sizeof( tree));
  new->a = a;
  new->b = b;
  new->parent = parent;

  new->left = left;
  new->right = right;
  new->hp = hp;
  new->t = NULL;
  new->gammas = malloc( 2*sizeof( double));
  new->gammas[0] = NAN;
  new->gammas[1] = NAN;
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

double hptree_get_ehp( tree *node, int r) {
  assert( node->hp);
  if( node->info.hp->lenehp < r || node->info.hp->ehp[r-1] == -1) {
    printf("Illegal access\n");
    exit(1);
  }
  return node->info.hp->ehp[r-1];
}

void hptree_set_ehp( tree *node, int r, double val) {
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

double hptree_get_tehp( tree *node, int r) {
  assert( node->hp);
  if( node->info.hp->lentehp < r || node->info.hp->tehp[r-1] == -1) {
    printf("Illegal access\n");
    exit(1);
  }
  return node->info.hp->tehp[r-1];
}

void hptree_set_tehp( tree *node, int r, double val) {
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

double hptree_get_e( tree *node, int r) {
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
    return (node->info.hp->e[r-1] = error( node, r));
  }
}

double htree_get_e( tree *node) {
  assert( !node->hp);
  if( node->info.h->e == -1) {
    node->info.h->e = error( node, node->info.h->r);
  }
  return node->info.h->e;
}

void tree_set_gamma( tree *node, int i, double gamma) {
  assert( node->lengammas > i);
  node->gammas[i] = gamma;
  if( i > node->hgammas) {
    node->hgammas = i;
  }
}

void tree_assert_gamma_len( tree *node, int i) {
  if( node->lengammas <= i) {
    int j;
    int pow2 = pow2roundup( i);
    node->gammas = realloc( node->gammas, 2*pow2*sizeof( double));
    for( j = node->hgammas+1; j < 2*pow2; j++) {
      node->gammas[j] = NAN;
    }
    node->lengammas = 2*pow2;
  }
}

double tree_get_gamma( tree *node, int i) {
  assert( node->lengammas > i);
  return node->gammas[i];
}

int tree_is_leaf( tree *node) {
  return node->left == NULL;
}

void tree_trim( workspace *w, tree *node) {
  assert( !tree_is_leaf( node));
  if( !tree_is_leaf( node->left)) {
    tree_trim( w, node->left);
  }
  workspace_remove_leaf( w, node->left);
  node->left = NULL;
  if( !tree_is_leaf( node->right)) {
    tree_trim( w, node->right);
  }
  workspace_remove_leaf( w, node->right);
  node->right = NULL;
  workspace_add_leaf( w, node);
}

void tree_subdivide( workspace *w, tree *node) {
  assert( tree_is_leaf( node));

  //point midway between element boundary
  double h = (node->a + node->b)/2;

  //create trees
  tree *lt = tree_create( node->a, h, node, NULL, NULL, node->hp);
  tree *rt = tree_create( h, node->b, node, NULL, NULL, node->hp);
  node->left = lt;
  node->right = rt;

  //change leaves
  workspace_remove_leaf( w, node);
  workspace_add_leaf( w, lt);
  workspace_add_leaf( w, rt);

  //set error
  if( !node->hp && node->info.h->r > 0) {
    lt->info.h->r = node->info.h->r;
    lt->info.h->te = inverror( htree_get_e( lt), node->info.h->te);
    rt->info.h->r = node->info.h->r;
    rt->info.h->te = inverror( htree_get_e( rt), node->info.h->te);
  }
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

/* Find inner nodes of the total partition. */
void tree_inner_nodes( workspace *w, tree ***inners, int *ninners) {
  tree **queue = malloc( w->nleaves * sizeof( tree *));
  *inners = malloc( w->nleaves * sizeof( tree *));
  *ninners = 0;
  int nqueue = 1;
  queue[0] = w->root;

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
