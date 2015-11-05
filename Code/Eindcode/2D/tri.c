#include <stdio.h>
#include <stdlib.h>
#include "tri.h"
#include "workspace.h"
#include "partition.h"
#include "triangulate.h"

/* definition 2.21 */
double tri_vol( workspace *w, tri *t) {
  point a = w->points[t->p[0]];
  point b = w->points[t->p[1]];
  point c = w->points[t->p[2]];
  return 0.5*( - b.x * a.y + c.x * a.y + a.x * b.y 
               - c.x * b.y - a.x * c.y + b.x * c.y);
}

/* Lemma 2.23 */
point tri_ref2t( workspace *w, tri *t, point p) {
  point a = w->points[t->p[0]];
  point b = w->points[t->p[1]];
  point c = w->points[t->p[2]];

  point Gp = {.x = p.x*(b.x - a.x) + p.y*(c.x - a.x) + a.x,
              .y = p.x*(b.y - a.y) + p.y*(c.y - a.y) + a.y};
  return Gp;
}

tri *tri_create( workspace *w, int a, int b, int c) {
  tri *t = malloc( sizeof( tri));
  t->p[0] = a;
  t->p[1] = b;
  t->p[2] = c;
  if( (t->vol = tri_vol( w, t)) < 0.0) {
    printf("wrong input: tri volume too small\n");
    exit(1);
  } else {
    return t;
  }
}

void tri_free( tri *t) {
  free( t);
}
