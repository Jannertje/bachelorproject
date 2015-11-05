#pragma once
#include "types.h"

workspace *workspace_init( tree *root);
void workspace_free( workspace *w);

/* various print functions */
void workspace_print( workspace *w);
void workspace_print_root( workspace *w);
void workspace_print_leaves( workspace *w);

/* print in a way suitable for plotting */
void workspace_print_plot( workspace *w);

/* mutate leaves of workspace */
int workspace_add_leaf( workspace *w, tree *node);
void workspace_remove_leaf( workspace *w, tree *node);
