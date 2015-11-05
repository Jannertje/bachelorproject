#pragma once
#include "types.h"

/* get node along edge */
tree *edge_get( tree **e, int j, int k);
/* reset edge value */
void edge_reset( tree **e, int j, int k);
/* set edge value */
void edge_set( tree **e, int j, int k, tree *node);

void edge_matrix_expand( workspace *w);
