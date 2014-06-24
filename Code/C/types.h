#pragma once
#include "main.h"

typedef struct point {
  double x, y;
} point;

typedef struct tri {
  int p[3]; //indices in points array; 0 is the newest vertex.
  double vol;
} tri;

typedef struct hp {
  int r;
  double te, q;

  double *e;
  int lene;

  double *ehp;
  int lenehp;

  double *tehp;
  int lentehp;

  double **coeffs;
  int lencoeffs;
} hp;

typedef struct h {
  int r;
  double e, te;
  double *coeffs;
} h;

typedef struct tree {
  int i; //index in tris array
  struct tree *parent, *left, *right;
  struct tree *t; //to overcome circular dependency

  union info {
    hp *hp;
    h *h;
  } info;
  int hp;
} tree;

typedef struct workspace {
  int npoly;
  int is_matched, is_conform;
  point *points;
  int npoints, lenpoints;

  tri **tris;
  int ntris, lentris;

  tree **roots;
  int nroots, lenroots;

  tree **leaves;
  int nleaves, lenleaves;

  tree **edges;
  int nedges, lenedges;
} workspace;


