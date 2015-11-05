#include "types.h"
#include "helper.h"

/* find smallest n such that 2^n >= x */
int pow2roundup( int x) {
  if (x < 0) {
    return 0;
  }
  --x;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return x+1;
}

/* See definition of n(r) in Section 2.3.4 */
int find_n( int p) {
  int q;
  for( q = 0;; q++) {
    if( (q+2)*(q+1)/2 > p) {
      return q-1;
    }
  }
}

/* See Remark 2.17 */
point F( point st) {
  point Fst = {.x = (1.0-st.y)*(1.0+st.x)/4.0, .y = (st.y + 1.0)/2.0};
  return Fst;
}
