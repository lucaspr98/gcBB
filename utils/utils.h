/**
  \file utils.h
  \brief Utility macros and functions.
  \internal Guilherme P. Telles, 2009-2017.
**/

#ifndef UTILSH
#define UTILSH

#include <zlib.h>

/**
   \brief A key/frequency pair.
**/
struct freq {
  unsigned k; ///<\brief key
  unsigned f; ///<\brief frequency
};

typedef struct freq freq;


/**
   \def   swap(a,b)
   \brief Swap.
*/
#define swap(a,b) do { typeof(a) aux_ = (a); (a) = (b); (b) = aux_; } while (0)


/**
   \def   min(a,b)
   \brief Minimum.
*/
#define min(a,b) ((a) < (b) ? (a) : (b))


/**
   \def   max(a,b)
   \brief Maximum.
*/
#define max(a,b) ((a) > (b) ? (a) : (b))


void die(char* format, ...);

int is_power2(unsigned x);
unsigned first_power2(unsigned x);
unsigned first_prime(unsigned m);

void mean_stddev(void* A, unsigned n, char type, double* mean, double* stddev);

struct timespec tsdiff(struct timespec start, struct timespec stop);

int cmpi(const void *a, const void *b);
int cmpu(const void *a, const void *b);
int cmps(const void *a, const void *b);
int cmpuc(const void *a, const void *b);
int cmp_freq_k(const void *a, const void *b);
int cmp_freq_fk(const void *a, const void *b);
int ccmp_freq_fk(const void *a, const void *b);

#endif
