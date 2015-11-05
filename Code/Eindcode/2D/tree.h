#pragma once
#include "types.h"

tree *tree_create( workspace *w, int tri_index, 
                   tree *parent, tree *left, tree *right, int hp);
void tree_free_shallow( tree *node); //free only the top node
void tree_free_deep( tree *node); //recursively free the whole tree

/* boolean value. 1 if node has no children, 0 else. */
int tree_is_leaf( tree *node);

void tree_subdivide( workspace *w, tree *node);
void tree_trim( workspace *w, tree *node);

/* sort **nodes list, biggest error tilde to smallest. See Section 1.1.4 */
void htree_sort_insertion( tree **nodes, int n);

/* functionality pertaining to the different (modified) errors of the node */
double hptree_get_ehp( workspace *w, tree *node, int r);
void hptree_set_ehp( workspace *w, tree *node, int r, double val);
double hptree_get_tehp( workspace *w, tree *node, int r);
void hptree_set_tehp( workspace *w, tree *node, int r, double val);
double hptree_get_e( workspace *w, tree *node, int r);
double htree_get_e( workspace *w, tree *node);

/* functionality pertaining to the coefficients of the polynomial 
 * of best approximation */
void tree_set_gamma( workspace *w, tree *node, int i, double gamma);
void tree_assert_gamma_len( workspace *w, tree *node, int i);
double tree_get_gamma( workspace *w, tree *node, int i);

/* is node on edge of domain? */
int tree_is_on_edge( workspace *w, tree *node);

/* does the node contain a hanging vertex? */
int tree_has_hanging_vertex( workspace *w, tree *node);
void tree_get_ancestors( workspace *w, tree *node, int *nancestors, tree ***ancestors);
