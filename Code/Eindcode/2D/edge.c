#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "edge.h"
#include "workspace.h"
#include "pair.h"

/* edge matrix.
 * This is a npoints x npoints matrix with at position (j,k) the node with edge 
 * with points (w->points[j], w->points[k]). This is _NOT_ symmetric.
 */

/* we store the edge matrix as an array and use studzik pairing function to 
 * find the array element. this is to be able to quickly add a row and a column 
 * to the array without too much reshuffling of the data structure.
 */


tree *edge_get( tree **e, int j, int k) {
  int pi;
  pair_szudzik( j, k, &pi);
  return e[pi];
}

void edge_reset( tree **e, int j, int k) {
  edge_set( e, j, k, NULL);
}

void edge_set( tree **e, int j, int k, tree *node) {
  int pi;
  pair_szudzik( j, k, &pi);
  e[pi] = node;
}

void edge_matrix_expand( workspace *w) {
  if( w->nedges == w->lenedges) {
    int l = w->lenedges*w->lenedges;
    w->edges = realloc( w->edges, 4*l*sizeof( tree *));

    //clean newly allocated memory
    int i;
    for( i = l; i < 4*l; i++) {
      w->edges[i] = NULL;
    }
    w->lenedges *= 2;
  }
  w->nedges += 1;
}
