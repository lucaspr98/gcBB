#ifndef UTILS_H
#define UTILS_H

#include <stdarg.h>
#include <errno.h>
#include <string.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

#ifndef UCHAR_SIZE
#define UCHAR_SIZE 256
#endif

#define END_MARKER '$'

#ifndef M64
        #define M64 0
#endif

#if M64
	typedef int64_t  int_t;
	typedef uint64_t uint_t;
	#define PRIdN	 PRId64
	#define U_MAX  UINT64_MAX
	#define I_MAX  INT64_MAX
	#define I_MIN  INT64_MIN
#else
	typedef int32_t  int_t;
	typedef uint32_t uint_t;
	#define PRIdN	 PRId32
	#define U_MAX  UINT32_MAX
	#define I_MAX	 INT32_MAX
	#define I_MIN	 INT32_MIN
#endif

typedef uint32_t int_text;

//Schemes for reversed string:
//1. T^r = a^r b^r c^r; 
//2. T^r = c^r b^r a^r
#define REVERSE_SCHEME 2 

/**********************************************************************/

#define swap(a,b) do { typeof(a) aux_a_b = (a); (a) = (b); (b) = aux_a_b; } while (0)

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

void   time_start(time_t *t_time, clock_t *c_clock);
double time_stop(time_t t_time, clock_t c_clock);

void die(const char* where);
void dies(const char* where, char* format, ...);

int_t print_int(int_t* A, int_t n);
int_t print_char(char* A, int_t n);
int_t min_range(int_t* A, int_t l, int_t r);

// used for DNA RLE encoding 
unsigned char map(unsigned char c);
unsigned char unmap(unsigned char c);
unsigned char rle(unsigned char c, unsigned char run);

int bwt(int_t sa, unsigned char* str);

/**********************************************************************/

int_t* cat_int(unsigned char** R, int k, int_t *n);
unsigned char* cat_char(unsigned char** R, int k, size_t *n);
unsigned char* cat_char_rev(unsigned char** R, int k, size_t *n);

void qsort2(void *array, size_t nitems, size_t size, int (*cmp)(void*,void*));

#endif
