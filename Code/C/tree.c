#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "macros.h"
#include "tree.h"

tree *tree_create( boundary b,
                   tree *left, tree *right, 
                   tree *parent,
                   error_info e) {

  assert( b.a < b.b); //[a,b] is a proper interval
  assert( left != NULL || right == NULL); //left == NULL => right == NULL

  tree *new = malloc( sizeof( tree));
  new->b = b;
  new->forest = malloc( 2 * sizeof( tree));
  new->parent = parent;
  new->e = e;

  new->forest[0] = left;
  new->forest[1] = right;
  
  return new;
}

void printer( tree *self, int indent) {
  if( self == NULL) return;
  for( int i = 0; i < indent; i++) {
    printf("  ");
  }
  printf("[%g,%g]", self->b.a, self->b.b);
  if( tree_is_leaf( self)) {
    printf(" (%g)\n", self->e.real_error);
  } else {
    printf("\n");
    printer( self->forest[0], indent + 1);
    printer( self->forest[1], indent + 1);
  }
}

void tree_print( tree *self) {
  printer( self, 0);
}

int tree_subdivide( tree *self) {
  assert( tree_is_leaf( self));

  double a = self->b.a;
  double b = self->b.b;
  double h = (a+b)/2.0;
  boundary left = {a, h};
  boundary right = {h, b};
  error_info e = {-1.0, -1.0};

  self->forest[0] = tree_create( left, NULL, NULL, self, e);
  self->forest[1] = tree_create( right, NULL, NULL, self, e);

  return 0;
}

tree_list *tree_leaves( tree *self) {
  if( tree_is_leaf( self)) {
    return tree_list_create( self);
  } else {
    tree_list *list0 = tree_leaves( self->forest[0]);
    tree_list *list1 = tree_leaves( self->forest[1]);
    tree_list *cur = list0;

    while( cur->next != NULL) {
      cur = cur->next;
    }

    cur->next = list1;
    return list0;
  }
}

tree_list *tree_inner_nodes( tree *self) {
  tree_list *list = NULL;
  tree_list *queue = tree_list_create( self);

  while( queue != NULL) {
    tree *node = tree_list_popleft( &queue);
    if( !tree_is_leaf( node)) {
      list = tree_list_insert_beginning( list, node);
      if( queue == NULL) {
        queue = tree_list_create( node->forest[0]);
      } else {
        tree_list_append( queue, node->forest[0]);
      }
      tree_list_append( queue, node->forest[1]);
    }
  }

  return list;
}

int tree_is_leaf( tree *self) {
  return self->forest[0] == NULL;
}

int tree_is_root( tree *self) {
  return self->parent == NULL;
}

int tree_height( tree *self) {
  if( tree_is_leaf( self)) {
    return 0;
  } else {
    return 1 + MAX( tree_height( self->forest[0]), 
                    tree_height( self->forest[1]));
  }
}

int tree_depth( tree *self) {
  if( self->parent != NULL) {
    return 1 + tree_depth( self->parent);
  } else {
    return 0;
  }
}

tree *tree_root( tree *self) {
  if( self->parent != NULL) {
    return tree_root( self->parent);
  } else {
    return self;
  }
}

tree **tree_siblings( tree *self) {
  return self->parent->forest;
}

void tree_free_subtree( tree *self) {
  if( !tree_is_leaf( self)) {
    tree_free_subtree( self->forest[0]);
    tree_free_subtree( self->forest[1]);
  }

  free( self->forest);
  free( self);
}

tree_list *tree_list_create( tree *node) {
  tree_list *list = malloc( sizeof( tree_list));
  list->node = node;
  list->next = NULL;
  return list;
}

tree_list *tree_list_insert_after( tree_list *list, tree *node) {
  tree_list *new = tree_list_create( node);
  new->next = list->next;
  list->next = new;
  return new;
}

tree_list *tree_list_append( tree_list *list, tree *node) {
  while( list->next != NULL) {
    list = list->next;
  }
  return tree_list_insert_after( list, node);
}

tree_list *tree_list_insert_beginning( tree_list *list, tree *node) {
  tree_list *new = tree_list_create( node);
  new->next = list;
  return new;
}

tree *tree_list_popleft( tree_list **list) {
  tree *node = (*list)->node;
  tree_list *tmp = *list;
  *list = (*list)->next;
  free( tmp);

  return node;
}

void tree_list_print( tree_list *list) {
  while( list != NULL) {
    printf("[%g,%g] ", list->node->b.a, list->node->b.b);
    list = list->next;
  }
  printf("\n");
}

int tree_list_remove( tree_list *list, tree_list *node) {
  while( list->next != NULL && list->next != node) {
    list = list->next;
  }

  if( list->next != NULL) {
    list->next = node->next;
    free( node);
    return 0;
  }
  return 1;
}

void tree_list_free( tree_list *list) {
  if( list == NULL) return;
  while( list->next != NULL) {
    tree_list_remove( list, list->next);
  }
  free( list);
}

