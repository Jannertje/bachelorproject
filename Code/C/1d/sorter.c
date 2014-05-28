#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <limits.h>

#include "tree.h"

int bins_next_pow2( double val) {
  return (int) log2( val);
}

int sorter_bins( tree_list *leaves, tree_list **bests) {
  int bin = INT_MIN; //highest bin (= power of 2) so far
  tree_list *cur = leaves;
  while( cur != NULL) {
    int curbin;
    error_info e = tree_error_info( cur->node);
    assert( e.error >= 0.0);
    curbin = bins_next_pow2( e.error);
    if( curbin > bin) {
      tree_list_free( *bests);
      *bests = tree_list_create( cur->node);
      bin = curbin;
    } else if( curbin == bin) {
      tree_list_append( *bests, cur->node);
    }

    cur = cur->next;
  }
  return 0;
}

/* TODO: sorteer je wel een lijst hier? het lijkt gewoon linear complexity te hebben!
 */
int sorter_sort( tree_list *leaves, tree_list **bests) {
  double best_error = -1.0;
  tree_list *cur = leaves;
  while( cur != NULL) {
    double error = tree_error_info( cur->node).error;
    assert( error >= 0.0);

    if( error > best_error) {
      tree_list_free( *bests);
      *bests = tree_list_create( cur->node);
      best_error = error;
    } else if( error == best_error) {
      tree_list_append( *bests, cur->node);
    }
    
    cur = cur->next;
  }
  return 0;
}
