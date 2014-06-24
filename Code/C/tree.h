#pragma once
#include "types.h"

tree *tree_create( workspace *w, int tri_index, tree *parent, tree *left, tree *right, int hp);
void tree_free_shallow( tree *node); //free only the top node
void tree_free_deep( tree *node); //recursively free the whole tree

int tree_is_leaf( tree *node);

void tree_subdivide( workspace *w, tree *node);
void tree_trim( workspace *w, tree *node);
void htree_sort_insertion( tree **nodes, int n);

double hptree_get_ehp( workspace *w, tree *node, int r);
void hptree_set_ehp( workspace *w, tree *node, int r, double val);
double hptree_get_tehp( workspace *w, tree *node, int r);
void hptree_set_tehp( workspace *w, tree *node, int r, double val);
double hptree_get_e( workspace *w, tree *node, int r);
double htree_get_e( workspace *w, tree *node);

void tree_set_gamma( workspace *w, tree *node, int i, double gamma);
double tree_get_gamma( workspace *w, tree *node, int i);

int tree_is_on_edge( workspace *w, tree *node);
int tree_has_hanging_vertex( workspace *w, tree *node);
