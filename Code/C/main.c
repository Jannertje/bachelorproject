#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "types.h"
#include "tri.h"
#include "tree.h"
#include "workspace.h"
#include "partition.h"
#include "triangulate.h"
#include "partition.h"
#include "treegen.h"
#include "main.h"
#include "pkd.h"

double f( double x, double y) {
  if( x+y < 1.1) {
    return 0;
  }
  else {
    return 1;
  }
}

int main( void) {
  int i;

  int hp = 1;

  workspace *w = workspace_init();

  point points[6] = {
    {.x = 0.5, .y = 1},
    {.x = 0.5, .y = 0.5},
    {.x = 0, .y = 0.5},
    {.x = 0, .y = 0}, 
    {.x = 1, .y = 0}, 
    {.x = 1, .y = 1},
  };
  w->npoly = sizeof( points)/sizeof( points[0]);
  for( i = 0; i < w->npoly; i++) {
    workspace_add_point( w, points[i]);
  }

  triangulate( w);
  partition_setup( w, hp);
  partition_match( w);

  if( hp) {
    treegen_hp( w, 7);
  } else {
    treegen_h( w, 10, 2);
  }
  partition_make_conform( w);
  //workspace_print_plot( w);
  for( i = 0; i < w->nleaves; i++) {
    int j;
    printf("%i\n", w->leaves[i]->i);
    for( j = 0; j < w->leaves[i]->info.hp->lencoeffs; j++) {
      int r;
      double *coeffs;
      if( !hptree_get_coeffs( w, w->leaves[i], j, &r, &coeffs)) {
        printf("\t%i ", r);
        int k;
        for( k = 0; k < r; k++) {
          printf("%g ", coeffs[k]);
        }
        printf("\n");
      }
    }
  }
  return 0;
  workspace_print( w);

  workspace_free( w);

  return 0;
}
