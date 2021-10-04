/*
  Guilherme P. Telles, 2009-2017.
*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <errno.h>

#include "utils.h"



/**
  \brief Print a message and terminate.

  This function receives an input in printf style, prints the string, prints
  strerror and exits with current errno (or 1 if errno is zero).

  \param format A printf-like format string and arguments.
**/
void die(char* format, ...) {

  int err = errno;
  va_list val;

  va_start(val,format);
  vprintf(format,val);
  va_end(val);

  if (err) {
    printf("%s.\n", strerror(err));
    exit(err);
  }
  exit(1);
}



/**
  \brief Check for a power of 2.
**/
int is_power2(unsigned x) {
  return ((x != 0) && !(x & (x - 1)));
}



/**
  \brief Get a power of 2 which is greater or equal to x.

  It will not work if unsigned is not 4 bytes long.

  \return Returns the first power of 2 that is greater than or equal to x.  On
  overflow it returns 0.
**/
unsigned first_power2(unsigned x) {

  // if (sizeof(unsigned) != 4) return 0;
  if (x == 0) return 1;

  // from Bit Twiddling Hacks http://graphics.stanford.edu/~seander/bithacks.html
  x--;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  x++;

  return x;
}



/**
  \brief Gets the smaller prime greater than or equal to m.

  \return Returns the smaller prime that is greater than or equal to m.  On
  overflow returns 0.
**/
unsigned first_prime(unsigned m) {

  if (m == 0 || m == 1) return 2;
  if (m == 2 || m == 3) return m;

  if (m%2 == 0)
    m += 1;

  while (m) {
    int i=0, s = (int) sqrt(m);

    if (m%3)
      for (i=6 ; i<=s; i+=6)
        if (m%(i+1) == 0 || m%(i-1) == 0)
          break;

    if (i > s+1)
      break;

    m+=2;
  }

  return m;
}



/**
  \brief Mean and the standard deviation of A[0..n-1].

  This is a straighforward two-pass implementation.

  \param type The type of elements in A.  Valid values are 'i' for int and 'd'
  for double.
**/
void mean_stddev(void* A, unsigned n, char type, double* mean, double* stddev) {

  int i;

  *mean = 0.0;

  if (type == 'i')
    for(i=0; i<n; i++)
      *mean += ((int*)A)[i];
  else if (type == 'd')
    for(i=0; i<n; i++)
      *mean += ((double*)A)[i];

  *mean /= ((double)n);

  *stddev = 0.0;

  if (type == 'i')
    for(i=0; i<n; i++)
      *stddev += (((int*)A)[i] - *mean)*(((int*)A)[i] - *mean);
  else if (type == 'd')
    for(i=0; i<n; i++)
      *stddev += (((double*)A)[i] - *mean)*(((double*)A)[i] - *mean);

  *stddev = sqrt(*stddev/((double)n));
}



/**
   \brief Evaluates the difference between two timespec structs.
**/
struct timespec tsdiff(struct timespec start, struct timespec stop) {

  struct timespec temp;

  temp.tv_sec = stop.tv_sec - start.tv_sec;
  temp.tv_nsec = stop.tv_nsec - start.tv_nsec;

  if (temp.tv_nsec < 0) {
    temp.tv_sec -= 1;
    temp.tv_nsec += 1e9;
  }

  return temp;
}


/**
   \brief Comparison function for unsigned char.
**/
int cmpuc(const void *a, const void *b) {
  return (*((unsigned char*)a) == *((unsigned char*)b) ? 0
          : (*((unsigned char*)a) > *((unsigned char*)b) ? 1 : -1));
}



/**
   \brief Comparison function for int.
**/
int cmpi(const void *a, const void *b) {
  return (*((int*)a) == *((int*)b) ? 0 : (*((int*)a) > *((int*)b) ? 1 : -1));
}



/**
   \brief Comparison function for unsigned.
**/
int cmpu(const void *a, const void *b) {
  return (*((unsigned*)a) == *((unsigned*)b) ? 0 : (*((unsigned*)a) > *((unsigned*)b) ? 1 : -1));
}



/**
   \brief Comparison function for null terminated string.
**/
int cmps(const void *a, const void *b) {
  return strcmp(*(const char**) a, *(const char**) b);
}



/**
   \brief Comparison function for freq by key.
**/
int cmp_freq_k(const void *a, const void *b) {
  return (((freq*)a)->k == ((freq*)b)->k ? 0 : (((freq*)a)->k > ((freq*)b)->k ? 1 : -1));
}



/**
   \brief Comparison function for freq by frequence and then by key.
**/
int cmp_freq_fk(const void *a, const void *b) {

  if (((freq*)a)->f == ((freq*)b)->f)
    return (((freq*)a)->k == ((freq*)b)->k ? 0 : (((freq*)a)->k > ((freq*)b)->k ? 1 : -1));
  else
    return (((freq*)a)->f > ((freq*)b)->f ? 1 : -1);
}



/**
   \brief Complement comparison function for freq by frequence and then by key.
**/
int ccmp_freq_fk(const void *a, const void *b) {

  if (((freq*)a)->f == ((freq*)b)->f)
    return (((freq*)a)->k == ((freq*)b)->k ? 0 : (((freq*)a)->k > ((freq*)b)->k ? -1 : 1));
  else
    return (((freq*)a)->f > ((freq*)b)->f ? -1 : 1);
}
