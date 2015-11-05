#include <math.h>
#include "pair.h"

/* see Section 0.2 for the definitions used. */

void pair_cantor( int x, int y, int *pi) {
  *pi = (x+y)*(x+y+1)/2 + y;
}

/* Because the inverse cantor pairing function took a lot 
 * of CPU cycles, cache the first 36 values */
void pair_invcantor( int *x, int *y, int pi) {
  switch( pi) {
    case 0: *x = 0; *y = 0; return;
    case 1: *x = 1; *y = 0; return;
    case 2: *x = 0; *y = 1; return;
    case 3: *x = 2; *y = 0; return;
    case 4: *x = 1; *y = 1; return;
    case 5: *x = 0; *y = 2; return;
    case 6: *x = 3; *y = 0; return;
    case 7: *x = 2; *y = 1; return;
    case 8: *x = 1; *y = 2; return;
    case 9: *x = 0; *y = 3; return;
    case 10: *x = 4; *y = 0; return;
    case 11: *x = 3; *y = 1; return;
    case 12: *x = 2; *y = 2; return;
    case 13: *x = 1; *y = 3; return;
    case 14: *x = 0; *y = 4; return;
    case 15: *x = 5; *y = 0; return;
    case 16: *x = 4; *y = 1; return;
    case 17: *x = 3; *y = 2; return;
    case 18: *x = 2; *y = 3; return;
    case 19: *x = 1; *y = 4; return;
    case 20: *x = 0; *y = 5; return;
    case 21: *x = 6; *y = 0; return;
    case 22: *x = 5; *y = 1; return;
    case 23: *x = 4; *y = 2; return;
    case 24: *x = 3; *y = 3; return;
    case 25: *x = 2; *y = 4; return;
    case 26: *x = 1; *y = 5; return;
    case 27: *x = 0; *y = 6; return;
    case 28: *x = 7; *y = 0; return;
    case 29: *x = 6; *y = 1; return;
    case 30: *x = 5; *y = 2; return;
    case 31: *x = 4; *y = 3; return;
    case 32: *x = 3; *y = 4; return;
    case 33: *x = 2; *y = 5; return;
    case 34: *x = 1; *y = 6; return;
    case 35: *x = 0; *y = 7; return;
  }
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
