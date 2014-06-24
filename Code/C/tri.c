#include <stdio.h>
#include <stdlib.h>
#include "tri.h"
#include "workspace.h"
#include "partition.h"
#include "triangulate.h"

double tri_vol( workspace *w, tri *t);

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

double tri_vol( workspace *w, tri *t) {
  point a = w->points[t->p[0]];
  point b = w->points[t->p[1]];
  point c = w->points[t->p[2]];
  return 0.5*( - b.x * a.y + c.x * a.y + a.x * b.y 
               - c.x * b.y - a.x * c.y + b.x * c.y);
}

point tri_t2ref( workspace *w, tri *t, point p) {
  double v = t->vol;
  point a = w->points[t->p[0]];
  point b = w->points[t->p[1]];
  point c = w->points[t->p[2]];
  point Gp = {.x = 0.5/v*((c.y - a.y)*(p.x - a.x) + (a.x - c.x)*(p.y - a.y)),
              .y = 0.5/v*((a.y - b.y)*(p.x - a.x) + (b.x - a.x)*(p.y - a.y))};
  return Gp;
}

point tri_ref2t( workspace *w, tri *t, point p) {
  point a = w->points[t->p[0]];
  point b = w->points[t->p[1]];
  point c = w->points[t->p[2]];

  point Gp = {.x = p.x*(b.x - a.x) + p.y*(c.x - a.x) + a.x,
              .y = p.x*(b.y - a.y) + p.y*(c.y - a.y) + a.y};
  return Gp;
}

/*
int main( void) {
  int i;
  workspace *w = workspace_init();

  point points[3] = {
    {.x = 0, .y = 0}, 
    {.x = 1, .y = 0}, 
    {.x = 0, .y = 1},
  };

  for( i = 0; i < sizeof( points)/sizeof( points[0]); i++) {
    workspace_add_point( w, points[i]);
  }
  workspace_add_tri( w, tri_create( w, 0, 1, 2));

  i = 0;
  tri *t = w->tris[i];
  point p = w->points[t->p[1]];
  //point p = {0, 0};
  point Gp = tri_t2ref( w, t, p);
  workspace_print( w);
  printf("%i: %g %g -> %g %g\n", i, p.x, p.y, Gp.x, Gp.y);

  return 0;
}
*/
