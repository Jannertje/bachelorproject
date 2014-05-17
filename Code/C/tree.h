#ifndef TREE_H_
#define TREE_H_

#include "types.h"

tree *tree_create( boundary b,
                   tree *left, tree *right, 
                   tree *parent,
                   error_info e);
void tree_free_subtree( tree *old);

int tree_subdivide( tree *self);

int tree_is_leaf( tree *self);
int tree_is_root( tree *self);
int tree_height( tree *self);
int tree_depth( tree *self);
tree *tree_root( tree *self);
tree **tree_siblings( tree *self);
void tree_print( tree *self);
void tree_list_print( tree_list *list);

tree_list *tree_leaves( tree *self);
tree_list *tree_inner_nodes( tree *self);

tree_list *tree_list_create( tree *node);
tree_list *tree_list_insert_after( tree_list *list, tree *node);
tree_list *tree_list_append( tree_list *list, tree *node);
tree_list *tree_list_insert_beginning( tree_list *list, tree *node);
tree *tree_list_popleft( tree_list **list);
int tree_list_remove( tree_list *list, tree_list *node);
void tree_list_free( tree_list *list);

#endif
