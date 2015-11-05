#pragma once
#include <types.h>

tri *tri_create( workspace *w, int a, int b, int c);
void tri_free( tri *node);

point tri_ref2t( workspace *w, tri *t, point p);
