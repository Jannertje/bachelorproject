#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
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
  //return sin(2*M_PI*x) * cos(2*M_PI*y);
  //return pow(x*x+y*y, -1/6.0);
  //return exp(-5 * (2*x-1)*(2*x-1)) * ((x-1)*(x-1)*x*x*(6*y*y+6*y-1)*(6*y*y+6*y-1) + 1600*y*y*(0.025 - 0.075*y + 0.05*y*y + x*x*(-0.75 + 2.25*y - 1.5*y*y) + x*(0.2 - 0.6*y + 0.4*y*y) + x*x*x*( 0.6 - 1.5*y + y*y))*(0.025 - 0.075 + 0.05*y*y + x*x*(-0.75 + 2.25*y - 1.5*y*y) + x*(0.2 - 0.6*y + 0.4*y*y) + x*x*x*( 0.6 - 1.5*y + y*y)));
  if( 2*y > x) {
    return 1;
  } else {
    return 0;
  }
}

int main( int argc, char **argv) {
  int N = 10;
  if( argc > 1) {
    N = atoi( argv[1]);
  }
  int i;

  int plot = 1;
  int hp = 1;
  workspace *w = workspace_init();

  /*
  point points[16] = {
    {0/3.0, 3/3.0}, {1/3.0, 3/3.0}, {2/3.0, 3/3.0}, {3/3.0, 3/3.0}, 
    {0/3.0, 2/3.0}, {1/3.0, 2/3.0}, {2/3.0, 2/3.0}, {3/3.0, 2/3.0}, 
    {0/3.0, 1/3.0}, {1/3.0, 1/3.0}, {2/3.0, 1/3.0}, {3/3.0, 1/3.0}, 
    {0/3.0, 0/3.0}, {1/3.0, 0/3.0}, {2/3.0, 0/3.0}, {3/3.0, 0/3.0}, 
  };
  */
  point points[6] = { {0, 0}, {0, -1}, {1, -1}, {1, 1}, {-1, 1}, {-1, 0}};
  //point points[3] = {{0, 0}, {1, 0}, {0, 1}};
  /*
  double r = 1;
  point points[5] = {
    {0, 0},
    {-r, 0},
    //{-r, -r},
    {0, -r},
    //{r, -r},
    {r, 0},
    //{r, r},
    {0, r}
  };
  */

  w->npoly = sizeof( points)/sizeof( points[0]);
  for( i = 0; i < w->npoly; i++) {
    workspace_add_point( w, points[i]);
  }

  triangulate( w);
  /*
  tri *tris[18] = {
    tri_create( w, 4, 5, 0), tri_create( w, 1, 0, 5),
    tri_create( w, 5, 6, 1), tri_create( w, 2, 1, 6),
    tri_create( w, 6, 7, 2), tri_create( w, 3, 2, 7),
    tri_create( w, 8, 9, 4), tri_create( w, 5, 4, 9),
    tri_create( w, 9, 10, 5), tri_create( w, 6, 5, 10),
    tri_create( w, 10, 11, 6), tri_create( w, 7, 6, 11),
    tri_create( w, 12, 13, 8), tri_create( w, 9, 8, 13),
    tri_create( w, 13, 14, 9), tri_create( w, 10, 9, 14),
    tri_create( w, 14, 15, 10), tri_create( w, 11, 10, 15),
  };
  */
  /*
  tri *tris[2] = { tri_create( w, 0, 1, 3), tri_create( w, 2, 3, 1)};
  for( i = 0; i < sizeof( tris)/sizeof( tris[0]); i++) {
    workspace_add_tri( w, tris[i]);
  }
  */
  partition_setup( w, hp);
  //partition_match( w);

  if( hp) {
    treegen_hp( w, N);
  } else {
    treegen_h( w, N, 1);
  }
  //partition_make_conform( w);
  if( plot) {
    workspace_print_plot( w);
  }
  //workspace_print( w);

  workspace_free( w);

  return 0;
}
