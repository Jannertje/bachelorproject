#pragma once
#include "main.h"

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
  double a, b; //boundary of element
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

  double norm_f; //\|f\|^2_{2,\node}. See Theorem 2.29
} tree;

typedef struct workspace {
  /* root of partition */
  tree *root;

  /* array of leaves in partition */
  tree **leaves;
  int nleaves, 
      lenleaves;
} workspace;
