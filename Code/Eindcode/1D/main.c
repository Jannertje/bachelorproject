#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "types.h"
#include "helper.h"
#include "tree.h"
#include "workspace.h"
#include "treegen.h"
#include "main.h"

double f( double x) {
  return pow( x, -0.3);
  if( x < 3) return exp(x)/exp(3);
  else return sin(M_PI *x/2);
}

int main( int argc, char **argv) {
  int N = 3;
  if( argc > 1) {
    N = atoi( argv[1]);
  }

  int plot = 1;
  int hp = 1;
  int r = 2;
  workspace *w = workspace_init( tree_create( 0, 5, NULL, NULL, NULL, hp));

  if( hp) {
    treegen_hp( w, N);
  } else {
    treegen_h( w, N, r);
  }
  //partition_make_conform( w);
  if( plot) {
    workspace_print_plot( w);
    if( hp) {
      print_total_error_hp( w);
    } else {
      print_total_error_h( w, r);
    }
  } else {
    workspace_print( w);
  }
  //workspace_print( w);

  workspace_free( w);

  return 0;
}
