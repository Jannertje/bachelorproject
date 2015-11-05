#include <stdio.h>
#include <math.h>
#include "tree.h"
#include "types.h"
#include "helper.h"

/* find smallest n such that 2^n >= x */
int pow2roundup( int x) {
  if (x < 0) {
    return 0;
  }
  --x;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return x+1;
}

void print_total_error_hp( workspace *w) {
  int dof = 0;
  double error = 0;
  int i;
  for( i = 0; i < w->nleaves; i++) {
    int r = w->leaves[i]->info.hp->r;
    dof += r;
    error += hptree_get_ehp( w->leaves[i], r);
  }
  fprintf( stderr, "%i\t%g\n", dof, error);
}

void print_total_error_h( workspace *w, int r) {
  int dof = 0;
  double error = 0;
  int i;
  for( i = 0; i < w->nleaves; i++) {
    dof += r;
    error += htree_get_e( w->leaves[i]);
  }
  fprintf( stderr, "%i\t%g\n", dof, error);
}
