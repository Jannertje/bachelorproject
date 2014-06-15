#pragma once
#include "types.h"

tree *tree_create( int tri_index, tree *parent, tree *left, tree *right);
void tree_free_shallow( tree *node); //free only the top node
void tree_free_deep( tree *node); //recursively free the whole tree

int tree_is_leaf( tree *node);

void tree_subdivide( workspace *w, tree *node);
