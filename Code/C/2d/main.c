#include <stdlib.h>
#include <stdio.h>
#include "types.h"
#include "tree.h"
#include "algo.h"
#include <math.h>

int main( void) {
  /*
  point a = { 0, 0};
  point b = { 1, 0};
  point c = { 1, 1};
  boundary bound = { a, b, c};
  error_info e = {-1.0, -1.0};

  tree *init = htree_create( bound, NULL, NULL, NULL, e);
  tree_subdivide( init);
  tree_subdivide( init->forest[0]);
  tree_subdivide( init->forest[0]->forest[0]);
  //algo_makeconform( init);
  tree_list_print( tree_leaves( init));
  printf("\n");
  */
  error_info e = {-1.0, -1.0};
  point zero = {0, 0};
  point b = {1.0, 0};
  point c = {1.0, 1};
  boundary t1 = { zero, b, c};
  tree *tree_t1 = htree_create( t1, NULL, NULL, NULL, e);
  tree_subdivide( tree_t1);
  tree_subdivide( tree_t1->forest[0]);
  tree_subdivide( tree_t1->forest[0]->forest[0]);
  tree_subdivide( tree_t1->forest[1]);
  tree_subdivide( tree_t1->forest[1]->forest[1]);
  tree_list_print( tree_leaves( tree_t1));
  printf("\n");
  return 0;
}
