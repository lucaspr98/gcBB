/**
  \file dynarray.h
  \brief A dynamic array.

  A dynamic array for int or void* that swells and shrinks by a factor of 2 as needed.

  push, pop, unshift and shift are amortized O(1).  get() is O(1) and put is O(n) when
  resizing is needed, but if i is kept within the array bounds, it is O(1).

  \internal g.p.telles, 2015,2017, sep 2019.
**/


#ifndef DYNARRAY_H
#define DYNARRAY_H


/**
   \brief The dynamic array.
**/
struct dynarray {

  char type;
  void* data;

  size_t first;
  size_t size;
  size_t cap;
  size_t min_cap;
};

typedef struct dynarray dynarray;


dynarray* da_alloc(char type, size_t capacity, char shrink);
void da_free(dynarray* A);

extern int da_size(dynarray* A);

int da_push(dynarray* A, ...);
int da_pop(dynarray* A, ...);
int da_unshift(dynarray* A, ...);
int da_shift(dynarray* A, ...);

int da_set(dynarray* A, size_t i, ...);
int da_get(dynarray* A, size_t i, ...);

#endif
