#include <stdio.h>
#include <math.h>
#include "types.h"
#include "macros.h"
#include "tree.h"
#include "algo.h"

int makeconform_collinear( point p, point a, point b) {
  //cross product small?
  return fabs((b.x - a.x) * (p.y - a.y) - (p.x - a.x) * (b.y - a.y)) < MINAREA;
}

tree_list *makeconform_find_neighbours( tree *self, tree *node) {
  tree_list *leaves = tree_leaves( self);
  tree_list *neighbours = NULL;
  tree_list *cur = leaves;
  while( cur != NULL) {
    //for each vertex: check if collinear with edge
    tree *curnode = cur->node;
    point list[9][3] = {
      { node->b.a, curnode->b.a, curnode->b.b},
      { node->b.a, curnode->b.b, curnode->b.c},
      { node->b.a, curnode->b.c, curnode->b.a},
      { node->b.b, curnode->b.a, curnode->b.b},
      { node->b.b, curnode->b.b, curnode->b.c},
      { node->b.b, curnode->b.c, curnode->b.a},
      { node->b.c, curnode->b.a, curnode->b.b},
      { node->b.c, curnode->b.b, curnode->b.c},
      { node->b.c, curnode->b.c, curnode->b.a},
    };
    int i;
    for( i = 0; i < 9; i++) {
      if( makeconform_collinear( list[i][0], list[i][1], list[i][2])) {
        //check if we are really neighbours
        if( makeconform_collinear( list[(i+3)%9][0], list[i][1], list[i][2])) {
          //(i:0,i+3:0) is collinear to (i:1,i:2)
          if( list[i][1].x == list[i][2].x) {
            //if( list[
          } else {

          }
        } else {
          //(i:0,i+6:0) is collinear to (i:1,i:2)
        }
        if( neighbours != NULL) {
          neighbours = tree_list_insert_beginning( neighbours, curnode);
        } else {
          neighbours = tree_list_create( curnode);
        }
        break;
      }
    }
    cur = cur->next;
  }

  tree_list_free( leaves);
  return neighbours;
}

int makeconform_within( double p, double q, double r) {
  return (q < p && p < r) || (q > p && p > r);
}
//is p on the line segment through ab?
int makeconform_is_point_on( point p, point a, point b) {
  if( a.x != b.x) {
    if( !makeconform_within( p.x, a.x, b.x)) return 0;
  } else {
    if( !makeconform_within( p.y, a.y, b.y)) return 0;
  }
  return makeconform_collinear( p, a, b);
}

//TODO: make quicker
int makeconform_ishanging( tree *self, tree *node) {
  tree_list *leaves = tree_leaves( self);
  tree_list *cur = leaves;
  int ret = 0;
  while( cur != NULL) {
    tree *curnode = cur->node;
    if( boundary_compare( curnode->b, node->b)) {
      if( makeconform_is_point_on( curnode->b.a, node->b.a, node->b.b) ||
          makeconform_is_point_on( curnode->b.a, node->b.b, node->b.c) ||
          makeconform_is_point_on( curnode->b.a, node->b.c, node->b.a)) {
        ret = 1;
        break;
      }

      if( makeconform_is_point_on( curnode->b.b, node->b.a, node->b.b) ||
          makeconform_is_point_on( curnode->b.b, node->b.b, node->b.c) ||
          makeconform_is_point_on( curnode->b.b, node->b.c, node->b.a)) {
        ret = 1;
        break;
      }

      if( makeconform_is_point_on( curnode->b.c, node->b.a, node->b.b) ||
          makeconform_is_point_on( curnode->b.c, node->b.b, node->b.c) ||
          makeconform_is_point_on( curnode->b.c, node->b.c, node->b.a)) {
        ret = 1;
        break;
      }
    }

    cur = cur->next;
  }
  tree_list_free( leaves);
  return ret;
}

//See NVBstevenson p. 249
int algo_makeconform( tree *self) {
  int bisections = 0;
  tree_list *leaves = tree_leaves( self);
  tree_list *cur = leaves;
  tree_list *M = NULL;
  while( cur != NULL) {
    if( makeconform_ishanging( self, cur->node)) {
      if( M != NULL) {
        M = tree_list_insert_beginning( M, cur->node);
      } else {
        M = tree_list_create( cur->node);
      }
    }
    cur = cur->next;
  }
  tree_list_free( leaves);

  tree *curtree = NULL;
  while( M != NULL) {
    curtree = tree_list_popleft( &M);
    tree_subdivide( curtree);
    //check if is hanging and put into M
    int i;
    for( i = 0; i < 2; i++) {
      if( makeconform_ishanging( self, curtree->forest[i])) {
        if( M != NULL) {
          M = tree_list_insert_beginning( M, curtree->forest[i]);
        } else {
          M = tree_list_create( curtree->forest[i]);
        }
      }
      tree_list *neighbours = makeconform_find_neighbours( self, curtree->forest[i]);
      tree_list *curn = neighbours;
      while( curn != NULL) {

        curn = curn->next;
      }
      tree_list_free( neighbours);
    }
    //did it create a hanging vertex in some neighbour of curtree?
  }
  return bisections;
}
