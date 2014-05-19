#ifndef TYPES_H_
#define TYPES_H_

typedef struct boundary {
  double a;
  double b;
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
typedef double ( *error)( function, boundary, int);
typedef int ( *sorter)( tree_list *, tree_list **);

typedef struct algo_info {
  function f;
  boundary b;
  sorter s;
  error e;
} algo_info;

#endif
