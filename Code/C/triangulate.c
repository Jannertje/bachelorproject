#include <stdlib.h>
#include <stdio.h>
#include "types.h"
#include "tri.h"
#include "workspace.h"
#include "partition.h"
#include "triangulate.h"

typedef struct point_list {
  point p;
  int pi;
  struct point_list *prev, *next;
} point_list;

typedef point vector;

double cross( vector v1, vector v2) {
  return v1.x * v2.y - v1.y * v2.x;
}

int sameside( point p1, point p2, point a, point b) {
  vector ba = {.x = b.x - a.x, .y = b.y - a.y};
  vector p1a = {.x = p1.x - a.x, .y = p1.y - a.y};
  vector p2a = {.x = p2.x - a.x, .y = p2.y - a.y};
  return (cross( ba, p1a) * cross( ba, p2a)) >= 0;
}

int point_in_tri( point p, point a, point b, point c) {
  return sameside( p, a, b, c) && sameside( p, b, a, c) && sameside( p, c, a, b);
}

//is the interior angle of angle abc < pi/2?
int is_reflex( point a, point b, point c) {
  vector ba = {.x = b.x - a.x, .y = b.y - a.y};
  vector cb = {.x = c.x - b.x, .y = c.y - b.y};
  return cross( ba, cb) >= 0;
}

void triangulate( workspace *w) {
  if( w->npoly < 3) {
    printf("not enough points to triangulate.\n");
    return;
  } else if( w->npoly == 3) {
    workspace_add_tri( w, tri_create( w, 0, 1, 2));
    return;
  }

  //set up doubly linked list
  point_list *first, *cur, *last;
  int i;
  for( i = 0; i < w->npoly; i++) {
    point_list *new = malloc( sizeof( point_list));
    new->prev = cur;
    new->p = w->points[i];
    new->pi = i;
    if( i == 0) {
      first = new;
    } else {
      cur->next = new;
    }
    cur = new;
  }
  cur->next = first;
  last = cur;
  first->prev = last;

  cur = first;
  while( cur != last) {
    if( is_reflex( cur->prev->p, cur->p, cur->next->p)) {
      point_list *curfirst = cur->next->next;
      point_list *curlast = cur->prev->prev;
      point_list *curcur = curfirst;
      int in_tri = 0;
      while( curcur != curlast) {
        if( point_in_tri( curcur->p, cur->prev->p, cur->p, cur->next->p)) {
          in_tri = 1;
          break;
        }
        curcur = curcur->next;
      }
      //curcur last
      if( !in_tri && !point_in_tri( curcur->p, cur->prev->p, cur->p, cur->next->p)) {
        //we have ear
        workspace_add_tri( w, tri_create( w, cur->prev->pi, cur->pi, cur->next->pi));
        cur->prev->next = cur->next;
        cur->next->prev = cur->prev;
        point_list *temp = cur;
        cur = cur->next;
        free( temp);
        continue;
      }
    }
    cur = cur->next;
  }
  //cur last
  if( is_reflex( cur->prev->p, cur->p, cur->next->p)) {
    point_list *curfirst = cur->next->next;
    point_list *curlast = cur->prev->prev;
    point_list *curcur = curfirst;
    int in_tri = 0;
    while( curcur != curlast) {
      if( point_in_tri( curcur->p, cur->prev->p, cur->p, cur->next->p)) {
        in_tri = 1;
        break;
      }
      curcur = curcur->next;
    }
    //curcur last
    if( !in_tri && !point_in_tri( curcur->p, cur->prev->p, cur->p, cur->next->p)) {
      //we have ear
      workspace_add_tri( w, tri_create( w, cur->prev->pi, cur->pi, cur->next->pi));
      cur->prev->next = cur->next;
      cur->next->prev = cur->prev;
      point_list *temp = cur;
      cur = cur->next;
      free( temp);
    }
  }
  if( cur->next->next->next == cur) { 
    //just 3 nodes left
    workspace_add_tri( w, tri_create( w, cur->prev->pi, cur->pi, cur->next->pi));
    free( cur->prev);
    free( cur->next);
    free( cur);
  }
}

/*
int main( void) {
  int i;
  workspace *w = workspace_init();

  point points[6] = {
    {.x = 1, .y = 0}, 
    {.x = 1, .y = 1},
    {.x = 2, .y = 1},
    {.x = 2, .y = 2},
    {.x = 0, .y = 2},
    {.x = 0, .y = 0}, 
  };

  for( i = 0; i < 6; i++) {
    workspace_add_point( w, points[i]);
  }

  triangulate( w);
  partition_setup( w);
  partition_match( w);
  //workspace_print_plot( w);
  workspace_free( w);

  //printf("%p %i\n", w->tris[0], w->ntris);
  return 0;
}
*/
