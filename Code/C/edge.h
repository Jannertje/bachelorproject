#pragma once
#include "types.h"

tree *edge_get( tree **e, int j, int k);
void edge_reset( tree **e, int j, int k);
void edge_set( tree **e, int j, int k, tree *node);

void edge_matrix_expand( workspace *w);
