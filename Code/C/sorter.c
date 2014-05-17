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
    assert( cur->node->e.error >= 0.0);
    int curbin = bins_next_pow2( cur->node->e.error);
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
    double error = cur->node->e.error;
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
