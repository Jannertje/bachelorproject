#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <assert.h>
#include <string.h>
#include "workspace.h"
#include "helper.h"
#include "tree.h"

#define CMPEPS 1E-16

workspace *workspace_init( tree *root) {
  workspace *w = malloc( sizeof( workspace));
  
  w->root = root;

  w->leaves = malloc( 2*sizeof( tree *));
  w->nleaves = 1;
  w->leaves[0] = root;
  w->lenleaves = 2;

  return w;
}

void workspace_free( workspace *w) {
  tree_free_deep( w->root);
  free( w->leaves);
  free( w);
}

int workspace_add_leaf( workspace *w, tree *node) {
  if( w->nleaves == w->lenleaves) {
    w->leaves = realloc( w->leaves, 2*w->lenleaves*sizeof( tree *));
    w->lenleaves *= 2;
  }
  w->leaves[w->nleaves] = node;
  return w->nleaves++;
}

void workspace_remove_leaf( workspace *w, tree *node) {
  int j;
  int i = -1;
  for( j = 0; j < w->nleaves; j++) {
    if( node->a == w->leaves[j]->a && node->b == w->leaves[j]->b) {
      i = j;
      break;
    }
  }
  assert( i > -1);
  if( i == --w->nleaves) {
    return;
  }

  for( j = i; j < w->nleaves; j++) {
    w->leaves[j] = w->leaves[j+1];
  }
}

void workspace_print_root( workspace *w) {
  //todo
}

void workspace_print_leaves( workspace *w) {
  int i;
  printf("leaves (n = %i, len = %i)\n", w->nleaves, w->lenleaves);
  for( i = 0; i < w->nleaves; i++) {
    printf("%i: a=%g b=%g", i, w->leaves[i]->a, w->leaves[i]->b);
    if( w->leaves[i]->hp) {
      printf(" r=%i, te=%g\n", w->leaves[i]->info.hp->r, w->leaves[i]->info.hp->te);
    } else {
      printf(" r=%i, te=%g, e=%g\n", w->leaves[i]->info.h->r, w->leaves[i]->info.h->te, w->leaves[i]->info.h->e);
    }
  }
  printf("\n");
}
void workspace_print( workspace *w) {
  printf("=================\n");
  workspace_print_root( w);
  workspace_print_leaves( w);
  printf("=================\n");
}

void workspace_print_plot( workspace *w) {
  int i,j;
  if( w->root->hp) {
    for( i = 0; i < w->nleaves; i++) {
      printf("%g %g %i", w->leaves[i]->a, w->leaves[i]->b, w->leaves[i]->info.hp->r);
      for( j = 0; j < w->leaves[i]->info.hp->r; j++) {
        printf(" %g", tree_get_gamma( w->leaves[i], j));
      }
      printf("\n");
    }
  } else {
    for( i = 0; i < w->nleaves; i++) {
      printf("%g %g %i", w->leaves[i]->a, w->leaves[i]->b, w->leaves[i]->info.h->r);
      for( j = 0; j < w->leaves[i]->info.h->r; j++) {
        printf(" %g", tree_get_gamma( w->leaves[i], j));
      }
      printf("\n");
    }
  }
  printf("\n");
}
