#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>

#include "tree.h"
#include "error.h"
#include "sorter.h"
#include "algo.h"

boundary b = {-1.0, 1.0};
double f( double x) {
  /*
  if( x < 0) {
    return x;
  } else {
    return x*x;
  }
  */
  return sin(x);
}

void run( void) {
  algo_info i = { &f, b, &sorter_sort, &error_2norm};
  tree *tree = algo_binev2007( &i, 4, 1);
  tree_free_subtree( tree);
}

void showUsage( char *argv0) {
  printf("Usage: %s [-s bins|hist] [-e 2norm] [-r #>0] [-n #] [-a greedy|binev2004|binev2007|binev2013]", argv0);
}

int main( int argc, char **argv) {
  gsl_set_error_handler_off();

  tree *( *algo)() = algo_binev2007;
  int r = 1, n = 4;
  algo_info i = { &f, b, &sorter_sort, &error_2norm};
  int is_hp = 0;

  int c;
  while( (c = getopt( argc, argv, "s:e:r:n:a:")) != -1) {
    switch( c) {
      case 's':
        if( strcmp( optarg, "bins") == 0) {
          i.s = &sorter_bins;
        } else if( strcmp( optarg, "hist") != 0) {
          printf("Unknown argument %s for option -s.\n", optarg);
          showUsage( argv[0]);
          return 1;
        }
        break;
      case 'e':
        if( strcmp( optarg, "2norm") != 0) {
          printf("Unknown argument %s for option -e.\n", optarg);
          showUsage( argv[0]);
          return 1;
        }
        break;
      case 'r':
        if( atoi( optarg) > 0) {
          r = atoi( optarg);
        } else {
          printf("Unknown argument %s for option -r.\n", optarg);
          showUsage( argv[0]);
          return 1;
        }
        break;
      case 'n':
        n = atoi( optarg);
        break;
      case 'a':
        if( strcmp( optarg, "greedy") == 0) {
          algo = algo_greedy;
        } else if( strcmp( optarg, "binev2004") == 0) {
          algo = algo_binev2004;
        } else if( strcmp( optarg, "binev2007") == 0) {
          algo = algo_binev2007;
        } else if( strcmp( optarg, "binev2013") == 0) {
          algo = algo_binev2013;
          is_hp = 1;
        } else {
          printf("Unknown argument %s for option -a.\n", optarg);
          showUsage( argv[0]);
          return 1;
        }
        break;
      default:
        printf("Unknown option %c.\n", optopt);
        showUsage( argv[0]);
        return 1;
    }
  }

  if( !is_hp) {
    tree *tree = algo( &i, r, n);
    tree_free_subtree( tree);
  } else {
    tree *tree = algo( &i, n);
    tree_free_subtree( tree);
  }
  return 0;
}
