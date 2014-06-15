#pragma once
#include "types.h"

int edge_get_leaf( int *e, int j, int k);
void edge_reset( int *e, int j, int k);
void edge_set( int *e, int j, int k, int i);

void edge_matrix_expand( workspace *w);
