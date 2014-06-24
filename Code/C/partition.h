#pragma once
#include "types.h"

void partition_setup( workspace *w, int hp);
void partition_match( workspace *w);
void partition_refine( workspace *w, tree **nodes, int n);
void partition_refine_node( workspace *w, tree *node);
void partition_make_conform( workspace *w);
void partition_inner_nodes( workspace *w, tree ***inners, int *ninners);
