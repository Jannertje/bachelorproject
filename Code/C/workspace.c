#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <assert.h>
#include <string.h>
#include "workspace.h"
#include "tree.h"
#include "edge.h"

#define CMPEPS 1E-7

workspace *workspace_init( void) {
  workspace *w = malloc( sizeof( workspace));

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
  w->edges = malloc( 4*sizeof( int));
  for( i = 0; i < 4; i++) {
    w->edges[i] = -1;
  }
  w->nedges = 0;
  w->lenedges = 2;

  return w;
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

void workspace_print_leaves( workspace *w) {
  int i;
  printf("leaves (n = %i, len = %i)\n", w->nleaves, w->lenleaves);
  for( i = 0; i < w->nleaves; i++) {
    printf("%i: tri index %i\n", i, w->leaves[i]->i);
  }
  printf("\n");
}

void workspace_print_edge_matrix( workspace *w) {
  int i, j;
  printf("edge matrix (n = %i, len = %i)\n", w->nedges, w->lenedges);
  for( i = 0; i < w->nedges * w->nedges; i++) {
    printf("%i ", w->edges[i]);
  }
  printf("\n");
  for( i = 0; i < w->nedges; i++) {
    for( j = 0; j < w->nedges; j++) {
      printf("%4i ", edge_get_leaf( w->edges, i, j));
    }
    printf("\n");
  }
  printf("\n");
}

void workspace_print( workspace *w) {
  printf("=================\n");
  workspace_print_points( w);
  workspace_print_tris( w);
  workspace_print_leaves( w);
  workspace_print_edge_matrix( w);
  printf("=================\n");
}
