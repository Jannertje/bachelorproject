#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <assert.h>
#include <string.h>
#include "workspace.h"
#include "tri.h"
#include "tree.h"
#include "edge.h"

#define CMPEPS 1E-7

workspace *workspace_init( void) {
  workspace *w = malloc( sizeof( workspace));

  w->npoly = 0;
  w->is_matched = 0;
  w->is_conform = 0;

  w->points = malloc( 2*sizeof( point));
  w->npoints = 0;
  w->lenpoints = 2;

  w->tris = malloc( 2*sizeof( tri *));
  w->ntris = 0;
  w->lentris = 2;
  
  w->roots = malloc( 2*sizeof( tree *));
  w->nroots = 0;
  w->lenroots = 2;

  w->leaves = malloc( 2*sizeof( tree *));
  w->nleaves = 0;
  w->lenleaves = 2;

  int i;
  w->edges = malloc( 4*sizeof( tree *));
  for( i = 0; i < 4; i++) {
    w->edges[i] = NULL;
  }
  w->nedges = 0;
  w->lenedges = 2;

  return w;
}

void workspace_free( workspace *w) {
  int i;
  free( w->points);

  for( i = 0; i < w->ntris; i++) {
    tri_free( w->tris[i]);
  }
  free( w->tris);

  for( i = 0; i < w->nroots; i++) {
    tree_free_deep( w->roots[i]);
  }
  free( w->roots);
  free( w->leaves);
  free( w->edges);
  free( w);
}

int workspace_add_point( workspace *w, point p) {
  //see if this point exists
  int i;
  for( i = w->npoints-1; i >= 0; i--) {
    if( !gsl_fcmp( p.x, w->points[i].x, CMPEPS) && 
        !gsl_fcmp( p.y, w->points[i].y, CMPEPS)) {
      return i;
    }
  }

  //else
  edge_matrix_expand( w);
  if( w->npoints == w->lenpoints) {
    w->points = realloc( w->points, 2*w->lenpoints*sizeof( point));
    w->lenpoints *= 2;
  }
  w->points[w->npoints] = p;
  return w->npoints++;
}

int workspace_add_tri( workspace *w, tri *t) {
  if( w->ntris == w->lentris) {
    w->tris = realloc( w->tris, 2*w->lentris*sizeof( tri *));
    w->lentris *= 2;
  }
  w->tris[w->ntris] = t;
  return w->ntris++;
}

int workspace_add_root( workspace *w, tree *node) {
  if( w->nroots == w->lenroots) {
    w->roots = realloc( w->roots, 2*w->lenroots*sizeof( tree *));
    w->lenroots *= 2;
  }
  w->roots[w->nroots] = node;
  return w->nroots++;
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
    if( node->i == w->leaves[j]->i) {
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

void workspace_print_points( workspace *w) {
  int i;
  printf("points (n = %i, len = %i)\n", w->npoints, w->lenpoints);
  for( i = 0; i < w->npoints; i++) {
    printf("%i: (%g,%g)\n", i, w->points[i].x, w->points[i].y);
  }
  printf("\n");
}

void workspace_print_tris( workspace *w) {
  int i;
  printf("tris (n = %i, len = %i)\n", w->ntris, w->lentris);
  for( i = 0; i < w->ntris; i++) {
    printf("%i: %i %i %i\n", i, w->tris[i]->p[0], w->tris[i]->p[1], w->tris[i]->p[2]);
  }
  printf("\n");
}

void workspace_print_roots( workspace *w) {
  int i;
  printf("roots (n = %i, len = %i)\n", w->nroots, w->lenroots);
  for( i = 0; i < w->nroots; i++) {
    printf("%i: tri=%i\n", i, w->roots[i]->i);
  }
  printf("\n");
}

void workspace_print_leaves( workspace *w) {
  int i;
  printf("leaves (n = %i, len = %i)\n", w->nleaves, w->lenleaves);
  for( i = 0; i < w->nleaves; i++) {
    printf("%i: tri=%i", i, w->leaves[i]->i);
    if( w->leaves[i]->hp) {
      printf(" r=%i, te=%g\n", w->leaves[i]->info.hp->r, w->leaves[i]->info.hp->te);
    } else {
      printf(" r=%i, te=%g, e=%g\n", w->leaves[i]->info.h->r, w->leaves[i]->info.h->te, w->leaves[i]->info.h->e);
    }
  }
  printf("\n");
}

void workspace_print_edge_matrix( workspace *w) {
  int i, j;
  printf("edge matrix (n = %i, len = %i)\n", w->nedges, w->lenedges);
  for( i = 0; i < w->nedges; i++) {
    for( j = 0; j < w->nedges; j++) {
      tree *leaf = edge_get( w->edges, i, j);
      if( leaf == NULL) {
        printf("%4i ", -1);
      } else {
        printf("%4i ", leaf->i);
      }
    }
    printf("\n");
  }
  printf("\n");
}

void workspace_print( workspace *w) {
  printf("=================\n");
  workspace_print_points( w);
  workspace_print_tris( w);
  workspace_print_roots( w);
  workspace_print_leaves( w);
  workspace_print_edge_matrix( w);
  printf("=================\n");
}

void workspace_print_plot( workspace *w) {
  int i, j;
  for( i = 0; i < w->npoints; i++) {
    printf("%g %g\n", w->points[i].x, w->points[i].y);
  }
  printf("\n");
  for( i = 0; i < w->nleaves; i++) {
    tri *t = w->tris[w->leaves[i]->i];
    printf("%i %i %i\n", t->p[0], t->p[1], t->p[2]);
    if( w->leaves[i]->hp) {
      printf("%i", w->leaves[i]->info.hp->r);
      int k;
      double *coeffs;
      hptree_get_coeffs( w, w->leaves[i], w->leaves[i]->info.hp->r, &k, &coeffs);
      for( j = 0; j < k; j++) {
        printf(" %g", coeffs[j]);
      }
    } else {
      printf("%i", w->leaves[i]->info.h->r);
      int k;
      double *coeffs;
      htree_get_coeffs( w, w->leaves[i], &k, &coeffs);
      for( j = 0; j < k; j++) {
        printf(" %g", coeffs[j]);
      }
    }
    printf("\n");
  }
  printf("\n");
}
