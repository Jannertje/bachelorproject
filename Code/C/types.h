#pragma once
typedef struct point {
  double x, y;
} point;

typedef struct tri {
  int p[3]; //indices in points array; 0 is the newest vertex.
  double vol;
} tri;

typedef struct tree {
  int i; //index in tris array
  struct tree *parent, *left, *right;
} tree;

typedef struct edge_matrix_element {
  int tri, leaf;
} edge_matrix_element;

typedef struct workspace {
  point *points;
  int npoints, lenpoints;

  tri **tris;
  int ntris, lentris;

  tree **roots;
  int nroots, lenroots;

  tree **leaves;
  int nleaves, lenleaves;

  int *edges;
  int nedges, lenedges;
} workspace;


