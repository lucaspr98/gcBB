/**
   \file arrays.h

   \brief Functions for arrays and ragged arrays.

   Throughout the functions in this file a lower triangular matrix (ltm) of
   order n is a ragged array with the following form:

   <pre>
   null
   M[1][0]
   M[2][0] M[2][1]
   M[3][0] M[3][1] M[3][2]
   ...
   M[n-1][0] M[n-1][1] M[n-1][2] ... M[n-1][n-2]
   </pre>

   As illustrated, this representation doesn't allocate the main diagonal
   elements and has a null pointer at row 0, preserving standard indexing.

   Throughout the functions in this file, the following letters are used for
   types of arrays and matrices.

   <pre>
   d double

   h: short
   i: int
   j: int8_t
   k: int16_t
   l: int32_t
   m: int64_t

   u: unsigned
   v: uint8_t
   w: uint16_t
   x: uint32_t
   y: uint64_t

   c: char
   s: string
   </pre>

   \internal Guilherme P. Telles, Jan 2010, 2017.
**/

#include "abbrevs.h"

#ifndef ARRAYSH
#define ARRAYSH


int a_load(void** A, char type, unsigned n, char* filename);
void a_print(void* A, char type, unsigned n);
void a_print_as_matrix(void* A, char type, unsigned n, unsigned m);

void a_reverse(void* A, char type, unsigned n);
int a_isequal(void* A, void* B, char type, unsigned n);


unsigned* a_ctou(char* S, unsigned n);
char* u8toa(u8 a, char* buffer);
char* u32toa(u32 a, char* buffer);
char* u64toa(u64 a, char* buffer);


void** m_alloc(char type, unsigned n, unsigned m);
void** m_allocz(char type, unsigned n, unsigned m);
void m_free(void** M, unsigned n);

int m_load(void*** M, char type, size_t n, size_t m, char* filename);
int m_load_sltm(void*** M, char type, size_t n, char* filename);

int m_write(void** A, char type, unsigned n, unsigned m, char* filename);
int m_writes(void** A, char type, unsigned n, unsigned m,
             char* column_sep, char* row_sep, char* filename);

void m_print(void** M, char type, size_t n, size_t m);
void m_prints(void** M, char type, size_t frow, size_t trow, size_t fcol, size_t tcol);


void** ltm_alloc(char type, unsigned n);
void ltm_free(void** M, unsigned n);

void** ltm_dup(void** S, char type, unsigned n);

int ltm_load(void*** M, char type, unsigned n, char* filename);
void ltm_print(void** A, char type, unsigned n);
int ltm_write(void** A, char type, unsigned n, char* filename);

#endif
