#include "types.h"
#include "tree.h"
#include "algo.h"

tree_list *algo_nvb( tree *tree, tree_list *marked) {
  tree_list *curmark = tree_list_copy( marked);
  tree_list_print( marked);
  tree_list_print( curmark);
  return curmark;
}
