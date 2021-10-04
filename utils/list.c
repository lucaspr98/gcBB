/*
  Guilherme P. Telles, 2009-2015.
*/

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <stdarg.h>

#include "list.h"



/**
   \brief Allocate an empty list.

   Allocate an empty list of given type.

   \param type Either 'c' (char), 'i' (int), 'l' (long), 'd' (double) or 'v' (void*).

   \return On success it returns a new list.  On failure it returns NULL and
   errno remains as set by malloc() on failure.
**/
list* list_alloc(char type) {

  list* L = (list*) malloc(sizeof(list));
  if (!L) return 0;

  // Sentinels:
  L->head = malloc(sizeof(list_node));
  if (!L->head) {
    free(L);
    return 0;
  }

  L->tail = malloc(sizeof(list_node));
  if (!L->tail) {
    free(L->head);
    free(L);
    return L;
  }

  L->head->prev = 0;
  L->head->next = L->tail;
  L->tail->prev = L->head;
  L->tail->next = 0;

  L->type = type;
  L->size = 0;

  L->pi = -1;

  return L;
}




/**
   \brief Release a list and its nodes.

   When data type is void*, data pointed to by nodes will not be freed.
**/
void list_free(list* L) {

  while (L->head) {
    L->tail = L->head;
    L->head = L->head->next;
    free(L->tail);
  }

  free(L);
}



/**
   \brief Empty a list.

   Release the list nodes.
**/
void list_clear(list* L) {

  while (L->head->next != L->tail) {
    list_node* p = L->head;
    L->head = L->head->next;
    free(p);
  }

  L->head->prev = 0;
  L->tail->prev = L->head;
  L->size = 0;
  L->p = 0;
  L->pi = -1;
}




/**
   \brief The number of nodes in the list.
**/
inline size_t list_size(list* L) {
  return L->size;
}




/**
   \brief Add n data items to the tail of L.

   \param n The number of data items to push.

   \param ... The data items.  Their type should match the type specified at
   list creation.  They will be pushed one by one, in the order they were given
   to the function.

   \return The number of nodes added to the list.
**/
int list_push(list* L, unsigned n, ...) {

  va_list va;
  va_start(va,n);
  int i;

  for (i=0; i<n; i++) {
    list_node* p = malloc(sizeof(list_node));
    if (!p) return i;

    switch (L->type) {
    case 'c': p->datum.c = (char) va_arg(va,int); break;
    case 'i': p->datum.i = va_arg(va,int); break;
    case 'l': p->datum.l = va_arg(va,long); break;
    case 'd': p->datum.d = va_arg(va,double); break;
    case 'v': p->datum.v = va_arg(va,void*);
    }

    p->prev = L->tail->prev;
    p->next = L->tail;
    L->tail->prev->next = p;
    L->tail->prev = p;
    L->size++;
  }

  va_end(va);
  return i;
}




/**
   \brief  Add n data items to the head of L.

   \param n The number of data items to inject.

   \param ... The data items.  Their type should match that specified at list
   creation.  They will be injected one by one in the order given to the
   function.

   \return The number of nodes added to the list.
**/
int list_inject(list* L, unsigned n, ...) {

  va_list va;
  va_start(va,n);
  int i;

  for (i=0; i<n; i++) {
    list_node* p = malloc(sizeof(list_node));
    if (!p) return i;

    switch (L->type) {
    case 'c': p->datum.c = (char) va_arg(va,int); break;
    case 'i': p->datum.i = va_arg(va,int); break;
    case 'l': p->datum.l = va_arg(va,long); break;
    case 'd': p->datum.d = va_arg(va,double); break;
    case 'v': p->datum.v = va_arg(va,void*);
    }

    p->prev = L->head;
    p->next = L->head->next;
    L->head->next->prev = p;
    L->head->next = p;
    L->size++;
  }

  if (L->pi >= 0)
    L->pi += i;

  va_end(va);
  return i;
}




/**
   \brief Remove n data items from the tail of a list.

   \param n The number of data items to pop.

   \param ... Addresses where the values in the removed nodes will be placed,
   which type should be pointer to the type specified at list creation.  List
   elements will be popped one by one and the returning addresses will be filled
   in the order given to the function.

   \return The number of nodes removed from the list.
**/
int list_pop(list* L, unsigned n, ...) {

  va_list va;
  va_start(va,n);
  int i;

  if (n > L->size) {
    n = L->size;
    L->pi = -1;
  }
  else
    if (L->pi > -1 && L->pi >= L->size-n)
      L->pi = -1;

  for (i=0; i<n; i++) {
    list_node *p = L->tail->prev;

    switch (L->type) {
    case 'c': *((char*) va_arg(va,int*)) = p->datum.c; break;
    case 'i': *(va_arg(va,int*)) = p->datum.i; break;
    case 'l': *(va_arg(va,long*)) = p->datum.l; break;
    case 'd': *(va_arg(va,double*)) = p->datum.d; break;
    case 'v': *(va_arg(va,void**)) = p->datum.v;
    }

    L->tail->prev = p->prev;
    L->tail->prev->next = L->tail;
    free(p);
    L->size--;
  }

  va_end(va);
  return i;
}




/**
   \brief Remove n data items from the head of a list.

   \param n The number of data items to eject.

   \param ... Addresses where the values in the removed nodes will be placed,
   which type should be a pointer to the type specified at list creation.  List
   elements will be ejected one by one and the returning addresses will be
   filled in the order given to the function.

   \return The number of nodes removed from the list.
**/
int list_eject(list* L, unsigned n, ...) {

  if (n > L->size) {
    n = L->size;
    L->pi = -1;
  }
  else {
    if (L->pi > -1 && L->pi < n)
      L->pi = -1;
    else
      L->pi -= n;
  }

  va_list va;
  va_start(va,n);

  int i;
  for (i=0; i<n; i++) {
    list_node *p = L->head->next;

    switch (L->type) {
    case 'c': *((char*) va_arg(va,int*)) = p->datum.c; break;
    case 'i': *(va_arg(va,int*)) = p->datum.i; break;
    case 'l': *(va_arg(va,long*)) = p->datum.l; break;
    case 'd': *(va_arg(va,double*)) = p->datum.d; break;
    case 'v': *(va_arg(va,void**)) = p->datum.v;
    }

    L->head->next = p->next;
    L->head->next->prev = L->head;
    free(p);
    L->size--;
  }

  va_end(va);
  return i;
}




/**
   \brief Insert n data items at positions pos,pos+1,...,pos+n-1 of a list.

   \param pos The position where the fisrt node will be inserted. The list head
   is at position 0.  If pos >= |L| then the data will be inserted after the
   list tail.

   \param n  The number of data items to insert.

   \param ...  The data items.  Their type should match the type specified at
   list creation.

   \return  The number of nodes added to the list.
**/
int list_insert(list* L, unsigned pos, unsigned n, ...) {

  if (pos > L->size)
    pos = L->size;

  // Make p point to the node at pos-1:
  int i;
  list_get(L,pos,1,&i);
  list_node* p = L->p->prev;

  va_list va;
  va_start(va,n);

  for (i=0; i<n; i++) {
    list_node* q = malloc(sizeof(list_node));
    if (!q) return i;

    switch (L->type) {
    case 'c': q->datum.c = (char) va_arg(va,int); break;
    case 'i': q->datum.i = va_arg(va,int); break;
    case 'l': q->datum.l = va_arg(va,long); break;
    case 'd': q->datum.d = va_arg(va,double); break;
    case 'v': q->datum.v = va_arg(va,void*);
    }

    q->prev = p;
    q->next = p->next;
    p->next->prev = q;
    p->next = q;
    p = q;

    L->pi++;
    L->size++;
  }

  va_end(va);
  return i;
}



/**
   \brief Set data items at pos, pos+1, ..., pos+n-1.

   A sequential list traversal through successive calls to this function will
   take linear time.

   \param pos The position of the first node to set. The list head is at
   position 0.  Only pre-existing positions will be set, up to the end of the
   list.  When data type is void*, pre-existing data pointed to by updated nodes
   will not be freed.

   \param n The number of data items to set.

   \param ... The data items, whose type should be pointer to the data type
   specified at list creation.

   \return The number of data items that were set.
**/
int list_set(list* L, unsigned pos, unsigned n, ...) {

  if (pos >= L->size) return 0;

  if (pos+n > L->size)
    n = L->size-pos;

  // Makes L->p point to the node at pos:
  if (L->pi <= pos)
    if (pos-L->pi <= L->size-pos-1) {
      //printf("from p %d\n",L->pi);
      for ( ; L->pi<pos; L->pi++)
        L->p = L->p->next;
    }
    else {
      L->p = L->tail->prev;
      //printf("from tail\n",L->pi);
      for (L->pi=L->size-1; L->pi>pos; L->pi--)
        L->p = L->p->prev;
    }
  else
    if (pos < L->pi-pos) {
      L->p = L->head->next;
      //printf("from head\n",L->pi);
      for (L->pi=0; L->pi<pos; L->pi++)
        L->p = L->p->next;
    }
    else {
      //printf("from p %d\n",L->pi);
      for ( ; L->pi>pos; L->pi--)
        L->p = L->p->prev;
    }

  va_list va;
  va_start(va,n);

  int i;
  for (i=0; i<n; i++) {

    switch (L->type) {
    case 'c': L->p->datum.c = (char) va_arg(va,int); break;
    case 'i': L->p->datum.i = va_arg(va,int); break;
    case 'l': L->p->datum.l = va_arg(va,long); break;
    case 'd': L->p->datum.d = va_arg(va,double); break;
    case 'v': L->p->datum.v = va_arg(va,void*);
    }

    if (i < n-1) {
      L->p = L->p->next;
      L->pi++;
    }
  }

  va_end(va);
  return i;
}




/**
   \brief Get data at pos,pos+1,...,pos+n-1 from a list.

   A sequential list traversal through successive calls to this function will
   take linear time.

   \param pos  The position of the first node to retrieve. The list head is at
   position 0.

   \param n  The number of data items to get.

   \param ...  The addresses where the data items will be placed, whose type
   should be pointer to the data type specified at list creation.

   \return  The number of data items retrieved from the list.
**/
int list_get(list* L, unsigned pos, unsigned n, ...) {

  if (pos+n > L->size)
    n = L->size-pos;

  // Makes L->p point to the node at pos:
  if (L->pi <= pos)
    if (pos-L->pi <= L->size-pos-1) {
      //printf("from p %d\n",L->pi);
      for ( ; L->pi<pos; L->pi++)
        L->p = L->p->next;
    }
    else {
      L->p = L->tail->prev;
      //printf("from tail\n",L->pi);
      for (L->pi=L->size-1; L->pi>pos; L->pi--)
        L->p = L->p->prev;
    }
  else
    if (pos < L->pi-pos) {
      L->p = L->head->next;
      //printf("from head\n",L->pi);
      for (L->pi=0; L->pi<pos; L->pi++)
        L->p = L->p->next;
    }
    else {
      //printf("from p %d\n",L->pi);
      for ( ; L->pi>pos; L->pi--)
        L->p = L->p->prev;
    }

  va_list va;
  va_start(va,n);

  int i;
  for (i=0; i<n; i++) {

    switch (L->type) {
    case 'c': *((char*) va_arg(va,int*)) = L->p->datum.c; break;
    case 'i': *(va_arg(va,int*)) = L->p->datum.i; break;
    case 'l': *(va_arg(va,long*)) = L->p->datum.l; break;
    case 'd': *(va_arg(va,double*)) = L->p->datum.d; break;
    case 'v': *(va_arg(va,void**)) = L->p->datum.v;
    }

    if (i < n-1) {
      L->p = L->p->next;
      L->pi++;
    }
  }

  va_end(va);
  return i;
}




/**
  \brief Create an array from a list.

  Create an array whose size equals the list size and has the same values in the
  list in the same order.

  \return On success returns the new array.  If the array could not be created
  then returns NULL and errno remains set by malloc() on failure.
**/
void* list_toa(list* L) {

  unsigned int i = 0;
  list_node* p = L->head->next;

  if (L->type == 'i') {
    int* A = malloc(L->size*sizeof(int));
    if (!A) return 0;

    while (p->next) {
      A[i++] = p->datum.i;
      p = p->next;
    }
    return (void*) A;
  }

  else if (L->type == 'l') {
    long* A = malloc(L->size*sizeof(long));
    if (!A) return 0;

    while (p->next) {
      A[i++] = p->datum.l;
      p = p->next;
    }
    return (void*) A;
  }

  else if (L->type == 'd') {
    double* A = malloc(L->size*sizeof(double));
    if (!A) return 0;

    while (p->next) {
      A[i++] = p->datum.d;
      p = p->next;
    }
    return (void*) A;
  }

  else if (L->type == 'c') {
    char* A = malloc(L->size*sizeof(char));
    if (!A) return 0;

    while (p->next) {
      A[i++] = p->datum.c;
      p = p->next;
    }

    return (void*) A;
  }

  else {  // L->type == 'v'
    void** A = malloc(L->size*sizeof(void*));
    if (!A) return 0;

    while (p->next) {
      A[i++] = p->datum.v;
      p = p->next;
    }

    return A;
  }

  return 0;
}




/**
   \brief Print list contents on a line.
**/
void list_print(list* L) {

  list_node* p = L->head->next;
  while (p->next) {
    if (L->type == 'i')
      printf("%d ",p->datum.i);
    else if (L->type == 'l')
      printf("%ld ",p->datum.l);
    else if (L->type == 'd')
      printf("%.4lf ",p->datum.d);
    else if (L->type == 'c')
      printf("%c ",p->datum.c);
    else  // L->type == 'v'
      printf("%p ",p->datum.v);

    p = p->next;
  }
  printf("\n");
}



/**
   \brief Print list contents in reverse order on a line.
**/
void list_printr(list* L) {

  list_node* p = L->tail->prev;
  while (p->prev) {
    if (L->type == 'i')
      printf("%d ",p->datum.i);
    else if (L->type == 'l')
      printf("%ld ",p->datum.l);
    else if (L->type == 'd')
      printf("%.4lf ",p->datum.d);
    else if (L->type == 'c')
      printf("%c ",p->datum.c);
    else  // L->type == 'v'
      printf("%p ",p->datum.v);

    p = p->prev;
  }
  printf("\n");
}
