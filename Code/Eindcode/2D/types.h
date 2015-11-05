#pragma once
#include "main.h"

typedef struct point {
  double x, y;
} point;

typedef struct tri {
  int p[3]; //indices in workspace points array; 0 is the newest vertex.
  double vol; //volume of tri. see Definition 2.21
} tri;

typedef struct hp {
  int r; //node-specific degree r(\node)
  double te, //\tilde e(\node)
         q; //q(\node). See definition 1.17

  double *e; //array \{ e_1(\node), \ldots, e_N(\node) \}
  int lene; //length of array

  double *ehp; //array \{ e_1^{hp}(\node), \ldots, e_N^{hp}(\node) \}
  int lenehp; //length of array

  double *tehp; //array \{ \tilde e_1^{hp}(\node), \ldots, \tilde e_N^{hp}(\node) \}
  int lentehp; //length of array
} hp;

typedef struct h {
  int r; //degree r
  double e, //e_r(\node)
         te; //\tilde e_r(\node)
} h;

typedef struct tree {
  int i; //index in tris array
  struct tree *parent, *left, *right;
  struct tree *t; //see definition 1.17

  union info {
    hp *hp;
    h *h;
  } info;
  int hp; //is it an h- or hp-tree?

  double *gammas; //array of \{ \gamma_{j,k}: 0 \leq j,k; j+k \leq N\}
  int hgammas, //highest gamma that is set
      lengammas; //allocated length of array

  int gen; //generation, i.e., length of path to root

  double norm_f; //\|f\|^2_{2,\node}. See Theorem 2.29
} tree;

typedef struct workspace {
  int npoly; //amount of points in polygon
  int is_matched, 
      is_conform; 

  /* points in partition */
  point *points;
  int npoints,
      lenpoints;

  /* triangles in partition */
  tri **tris;
  int ntris,
      lentris;

  /* roots in partition, i.e., nodes in initial partition */
  tree **roots;
  int nroots, 
      lenroots;

  /* array of leaves in partition */
  tree **leaves;
  int nleaves, 
      lenleaves;

  /* edge matrix */
  tree **edges;
  int nedges, 
      lenedges;
} workspace;


