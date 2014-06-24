#ifndef TYPES_H_
#define TYPES_H_

typedef struct point {
  double x, y;
} point;

typedef struct boundary {
  point a, b, c;
} boundary;

typedef struct error_info {
  double error;
  double real_error;
} error_info;

typedef struct htree {
  error_info e;
} htree;

typedef struct hptree {
  error_info e;
  int r;
} hptree;

typedef struct tree {
  boundary b;
  struct tree **forest;
  struct tree *parent;
  union t {
    htree *h;
    hptree *hp;
  } t;
  int is_hp;
} tree;

typedef struct tree_list {
  tree *node;
  struct tree_list *next;
} tree_list;

typedef double ( *function)( double x);
typedef double ( *error)( function, boundary, boundary, int);
typedef int ( *sorter)( tree_list *, tree_list **);

typedef struct algo_info {
  function f;
  boundary b;
  sorter s;
  error e;
} algo_info;

int point_compare( point p1, point p2);
int boundary_compare( boundary b1, boundary b2);
point point_halfway( point a, point b);
double boundary_area( boundary b);
void boundary_print( boundary b, int pretty);

#endif
