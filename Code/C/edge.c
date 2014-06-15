#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "edge.h"
#include "workspace.h"
#include "pair.h"

/* edge matrix.
 * This is a nleaves x nleaves matrix with at position (j,k) the index in leaves of the leaf with edge with points (w->points[j], w->points[k]). This is _NOT_ symmetric.
 */

/* we store the edge matrix as an array and use studzik pairing function to find the array element. this is to be able to quickly add a row and a column to the array without too much reshuffling of the data structure.
 */
int edge_get_leaf( int *e, int j, int k) {
  int pi;
  pair_szudzik( j, k, &pi);
  return e[pi];
}

void edge_reset( int *e, int j, int k) {
  edge_set( e, j, k, -1);
}

void edge_set( int *e, int j, int k, int i) {
  int pi;
  pair_szudzik( j, k, &pi);
  e[pi] = i;
}

void edge_matrix_expand( workspace *w) {
  if( w->nedges == w->lenedges) {
    int l = w->lenedges*w->lenedges;
    w->edges = realloc( w->edges, 4*l*sizeof( int));

    //clean newly allocated memory
    int i;
    for( i = l; i < 4*l; i++) {
      w->edges[i] = -1;
    }
    w->lenedges *= 2;
  }
  w->nedges += 1;
}
