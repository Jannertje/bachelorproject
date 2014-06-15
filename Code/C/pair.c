#include <math.h>
#include "pair.h"

void pair_cantor( int x, int y, int *pi) {
  *pi = (x+y)*(x+y+1)/2 + y;
}

void pair_invcantor( int *x, int *y, int pi) {
  int w = (int) floor((sqrt(8*pi+1)-1)/2);
  int t = (w*w + w)/2;
  *y = pi - t;
  *x = w - *y;
}

void pair_szudzik( int x, int y, int *pi) {
  if( x < y) {
    *pi = y*y+x;
  } else {
    *pi = x*x+x+y;
  }
}

void pair_invszudzik( int *x, int *y, int pi) {
  int w = (int) floor(sqrt(pi));
  int wsq = w*w;
  if( pi - wsq < w) {
    *x = pi - wsq;
    *y = w;
  } else {
    *x = w;
    *y = pi - wsq - w;
  }
}
