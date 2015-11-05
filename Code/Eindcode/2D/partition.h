#pragma once
#include "types.h"

/* given tris, create trees and setup initial partition */
void partition_setup( workspace *w, int hp);

/* create a matching partition by means of Algorithm A.11 */
void partition_match( workspace *w);

/* Create a conforming refinement where nodes[n] are refined. Algorithm A.13 */
void partition_refine( workspace *w, tree **nodes, int n);

/* Create a conforming refinement where node is refined. Algorithm A.5 */
void partition_refine_node( workspace *w, tree *node);

/* Algorithm A.7 */
void partition_make_conform( workspace *w);

/* Find inner nodes of the total partition. */
void partition_inner_nodes( workspace *w, tree ***inners, int *ninners);
