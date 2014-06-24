#pragma once
#include "types.h"

workspace *workspace_init( void);
void workspace_free( workspace *w);
void workspace_print( workspace *w);
void workspace_print_points( workspace *w);
void workspace_print_tris( workspace *w);
void workspace_print_roots( workspace *w);
void workspace_print_leaves( workspace *w);
void workspace_print_edge_matrix( workspace *w);

void workspace_print_plot( workspace *w);

int workspace_add_point( workspace *w, point p);
int workspace_add_tri( workspace *w, tri *t);
int workspace_add_root( workspace *w, tree *node);
int workspace_add_leaf( workspace *w, tree *node);
void workspace_remove_leaf( workspace *w, tree *node);
