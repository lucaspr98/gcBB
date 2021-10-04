/*
  Guilherme P. Telles, 2015,2017.
*/

#include <stdarg.h>
#include <stdlib.h>
#include <stdint.h>
#include <errno.h>
#include <string.h>

#include "dynarray.h"


/**
   \brief Allocate a dynamic array.

   \param type The type of data in the dynamic array.  Must be either 'i' (int)
   or 'v' (void*).

   \param capacity The initial capacity.  It is also the minimum capacity.

   \param shrink If true then the dynamic array will shrink.  If false it may
   grow but will never shrink.

   \return The address of a new dynamic array.  On failure it returns NULL and
   errno remains as set by malloc.
**/
dynarray* da_alloc(char type, size_t capacity, char shrink) {

  if (sizeof(char) != 1) {
    errno = EINVAL;
    return 0;
  }

  dynarray* A = malloc(sizeof(dynarray));
  if (!A) return 0;

  if (type == 'i')
    A->data = malloc(capacity*sizeof(int));
  else if (type == 'v')
    A->data = malloc(capacity*sizeof(void*));
  else
    return 0;

  if (!A->data) { free(A); return 0; }

  A->first = 0;
  A->cap = capacity;
  A->size = 0;
  A->type = type;

  if (shrink)
    A->min_cap = capacity;
  else
    A->min_cap = SIZE_MAX;

  return A;
}



/**
   \brief Release a dynamic array.
**/
void da_free(dynarray* A) {
  free(A->data);
  free(A);
}



/**
   \brief The number of items currently in the dynamic array.
**/
inline int da_size(dynarray* A) {
  return A->size;
}



/*
   Resize A to new_cap.
   Return 1 on success, 0 on malloc failure.
   The first array element will be at A->data[0] after resize.
*/
int da_resize(dynarray* A, size_t new_cap) {

  int s = 0;
  if (A->type == 'i')
    s = sizeof(int);
  else if (A->type == 'v')
    s = sizeof(void*);

  void* T = malloc(new_cap*s);
  if (!T) return 0;

  if (A->first+A->size >= A->cap) { // wraps:
    int n = A->cap - A->first;
    memcpy(T,((char*)A->data)+(A->first*s),n*s);
    memcpy(((char*)T)+(n*s),A->data,(A->size-n)*s);

    //if (A->type == 'i') {
    //  memcpy(T,((int*)A->data)+A->first,n*s);
    //  memcpy(((int*)T)+n,A->data,(A->size-n)*s);
    //}
    //else if (A->type == 'v') {
    //  memcpy(T,((void**)A->data)+A->first,n*s);
    //  memcpy(((void**)T)+n,((void**)A->data),(A->size-n)*s);
    //}
  }
  else // doesn't wrap:
    memcpy(T,((char*)A->data)+(A->first*s),A->size*s);

    //if (A->type == 'i')
    //  memcpy(T,((int*)A->data)+A->first,A->size*s);
    //else if (A->type == 'v')
    //  memcpy(T,((void**)A->data)+A->first,A->size*s);

  free(A->data);
  A->data = T;
  A->cap = new_cap;
  A->first = 0;

  return 1;
}



/**
   \brief Add a single data item to the end of the array.

   If the dynamic array is full, tries to double its capacity and then adds data.

   \param A A dynamic array.

   \param ... A single data item to be added, whose type is equal to the type
   specified to da_alloc().

   \return On success returns 1.  Whenever resizing the array is not possible it
   fails to push, returns 0 and leaves errno as set by malloc().
**/
int da_push(dynarray* A, ...) {

  if (A->size == A->cap)
    if (!da_resize(A,2*A->cap))
      return 0;

  va_list va;
  va_start(va,A);

  int p = A->first+A->size;
  if (p >= A->cap)
    p -= A->cap;

  if (A->type == 'i')
    *(((int*)A->data)+p) = va_arg(va,int);
  else if (A->type == 'v')
    *(((void**)A->data)+p) = va_arg(va,void*);

  va_end(va);
  A->size++;

  return 1;
}



/**
   \brief Remove a single data item from the end of the array.

   If the dynamic array is 1/4 full, this function halves its capacity and then
   removes data.  The capacity will never be smaller than the initial capacity.

   \param A A dynamic array.

   \param ... The address where the removed value will be placed, whose type is
   pointer to the type specified to da_alloc().

   \return On success it returns 1.  If the dynamic array is empty it returns 0.
**/
int da_pop(dynarray* A, ...) {

  if (!A->size) return 0;

  if (A->cap > A->min_cap && A->size <= A->cap/4)
    da_resize(A,A->cap/2);

  va_list va;
  va_start(va,A);

  A->size--;

  int p = A->first+A->size;
  if (p >= A->cap)
    p -= A->cap;

  if (A->type == 'i')
    *(va_arg(va,int*)) = *(((int*)A->data)+p);
  else if (A->type == 'v')
    *(va_arg(va,void**)) = *(((void**)A->data)+p);

  va_end(va);
  return 1;
}



/**
   \brief \brief Add a single data item to the beginning of the array.

   If the dynamic array is full, it tries to double its capacity and then adds
   data.

   \param A The dynamic array.

   \param ... A single data item to be added, whose type should be pointer to
   the type specified at dynamic array creation.

   \return On success returns 1.  When resizing the array is not possible it
   fails to unshift data, returns 0 and leaves errno as set by malloc().
**/
int da_unshift(dynarray* A, ...) {

  if (A->size == A->cap)
    if (!da_resize(A,2*A->cap))
      return 0;

  if (A->first == 0)
    A->first = A->cap-1;
  else
    A->first--;

  va_list va;
  va_start(va,A);

  if (A->type == 'i')
    *(((int*)A->data)+A->first) = va_arg(va,int);
  else if (A->type == 'v')
    *(((void**)A->data)+A->first) = va_arg(va,void*);

  va_end(va);

  A->size++;
  return 1;
}



/**
   \brief Remove a data item from the beginning of the array.

   If the dynamic array is 1/4 full, it halves its capacity and then removes data.
   The capacity will never be smaller than the initial capacity.

   \param A The dynamic array.

   \param ... The address where the removed value will be placed, whose type is
   pointer to the type specified to da_alloc().

   \return On success it returns 1.  If the dynamic array is empty it returns 0.
**/
int da_shift(dynarray* A, ...) {

  if (!A->size) return 0;

  if (A->cap > A->min_cap && A->size == A->cap/4)
    da_resize(A,A->cap/2);

  va_list va;
  va_start(va,A);

  if (A->type == 'i')
    *(va_arg(va,int*)) = *(((int*)A->data)+A->first);
  else if (A->type == 'v')
    *(va_arg(va,void**)) = *(((void**)A->data)+A->first);

  va_end(va);

  A->first++;
  if (A->first == A->cap)
    A->first = 0;

  A->size--;
  return 1;
}



/**
   \brief Set A[i].

   If i is outside A, tries to swell the array first.

   \param A A dynamic array.
   \param i The index.
   \param ... A single data item to set, whose type is equal to the type
   specified to da_alloc().

   \return On success returns 1.  Whenever resizing the array is not possible it
   leaves errno as set by malloc(), A remains unchanged and it returns 0.
**/
int da_set(dynarray* A, size_t i, ...) {

  if (i >= A->cap) {
    size_t k = A->cap;
    while (i >= k)
      k *= 2;

    if (!da_resize(A,k))
      return 0;
  }

  if (i >= A->size)
    A->size = i+1;

  size_t p = (A->first+i);
  if (p >= A->cap)
    p -= A->cap;

  va_list va;
  va_start(va,i);

  if (A->type == 'i')
    *(((int*)A->data)+p) = va_arg(va,int);
  else if (A->type == 'v')
    *(((void**)A->data)+p) = va_arg(va,void*);

  va_end(va);
  return 1;
}



/**
   \brief Get A[i].

   \param A A dynamic array.
   \param i The index.
   \param ... The address where A[i] will be copied to.  The type should be
   pointer to the type specified to da_alloc().

   \return On success returns 1.  If i is outside A it returns 0.
**/
int da_get(dynarray* A, size_t i, ...) {

  if (i >= A->size)
    return 0;

  size_t p = (A->first+i);
  if (p >= A->cap)
    p -= A->cap;

  va_list va;
  va_start(va,i);

  if (A->type == 'i')
    *(va_arg(va,int*)) = *(((int*)A->data)+p);
  else if (A->type == 'v')
    *(va_arg(va,void**)) = *(((void**)A->data)+p);

  va_end(va);
  return 1;
}
