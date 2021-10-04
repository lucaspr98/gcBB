/**
  \file list.h
  \brief A list.

  A list that may store either char, int, long, double or void*.

  This is a doubly-linked straightforward implementation, minimally tuned for
  sequential traversals.

  \internal Guilherme P. Telles, Oct 2009-2020.
**/

#ifndef LISTH
#define LISTH


typedef struct list_node {

  union {
    char c;
    int i;
    long l;
    double d;
    void* v;
  } datum;

  struct list_node* next;
  struct list_node* prev;

} list_node;




/**
   \brief The list.
**/
struct list {

  char type;
  list_node* head;
  list_node* tail;
  unsigned size;

  list_node* p;  // the most recently accessed node and
  int pi;        // its rank.

};

typedef struct list list;


list* list_alloc(char type);
void list_free(list* L);
void list_clear(list* L);

extern size_t list_size(list* L);

int list_push(list* L, unsigned n, ...);
int list_pop(list* L, unsigned n, ...);

int list_inject(list* L, unsigned n, ...);
int list_eject(list* L, unsigned n, ...);

int list_insert(list* L, unsigned pos, unsigned n, ...);

int list_get(list* L, unsigned pos, unsigned n, ...);
int list_set(list* L, unsigned pos, unsigned n, ...);

void* list_toa(list* L);

void list_print(list* L);
void list_printr(list* L);

#endif
