#ifndef TREE_H_
#define TREE_H_

#include "types.h"

int location_compare( location l1, location l2);

tree *htree_create( location l,
                   tree *left, tree *right, 
                   tree *parent,
                   error_info e);
error_info htree_error_info( tree *self);
int htree_r_node( tree *node, tree *TN, tree **node_in_TN);
tree *htree_copy( tree *self);

tree *hptree_create( location l,
                   tree *left, tree *right, 
                   tree *parent,
                   error_info e,
                   int r);
error_info hptree_error_info( tree *self);
int hptree_r( tree *self);


int tree_subdivide( tree *self);
void tree_trim_subtree( tree *old);
void tree_free_subtree( tree *old);

int tree_is_leaf( tree *self);
int tree_is_root( tree *self);
int tree_height( tree *self);
int tree_depth( tree *self);
tree *tree_root( tree *self);
tree **tree_siblings( tree *self);
void tree_print( tree *self);
void tree_list_print( tree_list *list);
int tree_num_leaves( tree *self);
tree *tree_find_node( tree *node, tree *tree);

error_info tree_error_info( tree *self);

tree_list *tree_nodes( tree *self);
tree_list *tree_leaves( tree *self);
tree_list *tree_inner_nodes( tree *self);

tree_list *tree_list_create( tree *node);
tree_list *tree_list_insert_after( tree_list *list, tree *node);
tree_list *tree_list_append( tree_list *list, tree *node);
tree_list *tree_list_insert_beginning( tree_list *list, tree *node);
tree *tree_list_popleft( tree_list **list);
tree *tree_list_find_node( tree *node, tree_list *list);
int tree_list_remove( tree_list *list, tree_list *node);
void tree_list_free( tree_list *list);

#endif
