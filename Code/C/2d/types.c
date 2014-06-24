#include <math.h>
#include <stdio.h>
#include "types.h"

int point_compare( point p1, point p2) {
  if( p1.x < p2.x) return -1;
  if( p1.x > p2.x) return 1;
  if( p1.y < p2.y) return -1;
  if( p1.y > p2.y) return 1;
  return 0;
}

/*
 * Returns a value in {-1,0,1} in accorance with strcmp
 */
int boundary_compare( boundary b1, boundary b2) {
  return point_compare( b1.a, b2.a) || point_compare( b1.b, b2.b) || point_compare( b1.c, b2.c);
}

point point_halfway( point a, point b) {
  point p = {.x = (a.x + b.x)/2, .y = (a.y + b.y)/2};
  return p;
}

double boundary_area( boundary bound) {
  point a = bound.a;
  point b = bound.b;
  point c = bound.c;
  double total = 0.5 * (-b.x*a.y + c.x*a.y + a.x*b.y - c.x*b.y - a.x*c.y + b.x*c.y); //signed triangle area http://mathworld.wolfram.com/TriangleArea.html
  return total;
}

void boundary_print( boundary b, int pretty) {
  if( pretty) {
    printf("(%g,%g), (%g,%g), (%g,%g)", b.a.x, b.a.y, b.b.x, b.b.y, b.c.x, b.c.y);
  } else {
    printf("%g %g %g %g %g %g", b.a.x, b.a.y, b.b.x, b.b.y, b.c.x, b.c.y);
  }
}
