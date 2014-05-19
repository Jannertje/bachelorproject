#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "macros.h"
#include "tree.h"

/*
 * Returns a value in {-1,0,1} in accorance with strcmp
 */
int boundary_compare( boundary b1, boundary b2) {
  if( b1.a < b2.a) return -1;
  if( b1.a > b2.a) return 1;
  if( b1.b < b2.b) return -1;
  if( b1.b > b2.b) return 1;
  return 0;
}

tree *copier( tree *self, tree *parent) {
  if( self == NULL) return NULL;
  assert( !self->is_hp);
  tree *new;
  new = htree_create( self->b, NULL, NULL, parent, self->t.h->e);
  new->forest[0] = copier( self->forest[0], new);
  new->forest[1] = copier( self->forest[1], new);
  return new;
}

tree *htree_copy( tree *self) {
  return copier( self, NULL);
}

tree *tree_create( boundary b,
                   tree *left, tree *right, 
                   tree *parent) {

  assert( b.a < b.b); //[a,b] is a proper interval
  assert( left != NULL || right == NULL); //left == NULL => right == NULL

  tree *new = malloc( sizeof( tree));
  new->b = b;
  new->forest = malloc( 2 * sizeof( tree));
  new->parent = parent;

  new->forest[0] = left;
  new->forest[1] = right;
  
  return new;
}

tree *hptree_create( boundary b,
                   tree *left, tree *right, 
                   tree *parent,
                   error_info e,
                   int r) {
  
  tree *new = tree_create( b, left, right, parent);

  new->is_hp = 1;
  new->t.hp = malloc( sizeof( hptree));
  new->t.hp->e = e;
  return new;
}

tree *htree_create( boundary b,
                   tree *left, tree *right, 
                   tree *parent,
                   error_info e) {
  
  tree *new = tree_create( b, left, right, parent);

  new->is_hp = 0;
  new->t.h = malloc( sizeof( htree));
  new->t.h->e = e;
  return new;
}

void printer( tree *self, int indent) {
  if( self == NULL) return;
  int i;
  for( i = 0; i < indent; i++) {
    printf("  ");
  }
  printf("%p:[%g,%g]", self, self->b.a, self->b.b);
  if( tree_is_leaf( self)) {
    if( !self->is_hp) {
      printf(" (%g)\n", htree_error_info( self).real_error);
    } else {
      printf(" (%g, %i)\n", hptree_error_info( self).real_error, hptree_r( self));
    }
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

  if( !self->is_hp) {
    self->forest[0] = htree_create( left, NULL, NULL, self, e);
    self->forest[1] = htree_create( right, NULL, NULL, self, e);
  }

  return 0;
}

int tree_num_leaves( tree *self) {
  if( tree_is_leaf( self)) {
    return 1;
  } else {
    int num0 = tree_num_leaves( self->forest[0]);
    int num1 = tree_num_leaves( self->forest[1]);
    return num0 + num1;
  }
}

/*
 * TODO: of course there are multiple traversals possible. choose the best for the job?
 */
tree_list *tree_nodes( tree *self) {
  tree_list *list = tree_list_create( self);
  if( tree_is_leaf( self)) {
    return list;
  } else {
    tree_list *list0 = tree_nodes( self->forest[0]);
    tree_list *list1 = tree_nodes( self->forest[0]);
    tree_list *cur = list0;
    while( cur->next != NULL) {
      cur = cur->next;
    }
    cur->next = list1;
    list->next = list0;
    return list;
  }
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

error_info htree_error_info( tree *self) {
  assert( !self->is_hp);
  return self->t.h->e;
}

error_info hptree_error_info( tree *self) {
  assert( self->is_hp);
  return self->t.hp->e;
}

error_info tree_error_info( tree *self) {
  if( self->is_hp) {
    return hptree_error_info( self);
  } else {
    return htree_error_info( self);
  }
}

int hptree_r( tree *self) {
  assert( self->is_hp);
  return self->t.hp->r;
}

tree *tree_find_node( tree *node, tree *self) {
  //pre-order traverse for the moment. TODO
  if( boundary_compare( self->b, node->b) == 0) {
    return self;
  }
  if( tree_is_leaf( self)) {
    return NULL;
  }
  tree *left = tree_find_node( node, self->forest[0]);
  if( left != NULL) {
    return left;
  }
  return tree_find_node( node, self->forest[1]);
}

int htree_r_node( tree *node, tree *TN, tree **node_in_TN) {
  assert( !node->is_hp);
  *node_in_TN = tree_find_node( node, TN);
  assert( *node_in_TN != NULL);
  return tree_num_leaves( *node_in_TN);
}

void tree_trim_subtree( tree *self) {
  if( !tree_is_leaf( self)) {
    tree_free_subtree( self->forest[0]);
    tree_free_subtree( self->forest[1]);
  }

  self->forest[0] = NULL;
  self->forest[1] = NULL;
}

void tree_free_subtree( tree *self) {
  if( !tree_is_leaf( self)) {
    tree_free_subtree( self->forest[0]);
    tree_free_subtree( self->forest[1]);
  }

  if( self->is_hp) {
    free( self->t.hp);
  } else {
    free( self->t.h);
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

tree *tree_list_find_node( tree *node, tree_list *list) {
  while( list != NULL) {
    if( boundary_compare( list->node->b, node->b) == 0) {
      return list->node;
    }
    list = list->next;
  }
  return NULL;
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

