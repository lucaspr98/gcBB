/*
  Guilherme P. Telles, 2010-2020.
*/

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#include "abbrevs.h"
#include "arrays.h"



/**
   \brief Reverse an array.

    a_reverse() reverses the n leading elements of an array of given type.

   \param type The type of *A.  Valid types are 'i' (int), or 'd' (double).
**/
void a_reverse(void* A, char type, unsigned n) {

  unsigned i;

  if (type == 'i') {
    int* a = (int*) A;
    for (i=0; i<n/2; i++) {
      int aux = a[i];
      a[i] = a[n-i-1];
      a[n-i-1] = aux;
    }
  }
  else if (type == 'd') {
    double* a = (double*) A;
    for (i=0; i<n/2; i++) {
      int aux = a[i];
      a[i] = a[n-i-1];
      a[n-i-1] = aux;
    }
  }
}




/**
   \brief Test arrays for equality.

   a_isequal() tests if the n leading elements of two arrays are equal.

   \param type The type of *A.  Valid types are 'i' (int) or 'd' (double).
   \return If A[0,n-1] = B[0,n-1] it returns 1, otherwise 0.
**/
int a_isequal(void* A, void* B, char type, unsigned n) {

  unsigned i;

  if (type == 'i') {
    int* a = (int*) A;
    int* b = (int*) B;
    for (i=0; i<n; i++)
      if (a[i] != b[i])
        return 0;
  }
  else if (type == 'd') {
    double* a = (double*) A;
    double* b = (double*) B;
    for (i=0; i<n; i++)
      if (a[i] != b[i])
        return 0;
  }

  return 1;
}




/**
   \brief Build an unsigned array from a char array.

   a_ctou() builds an unsigned array from the n leading chars in a char array.

   \returns The address of a new unsigned array with n elements whose values are
   equal to the first n values in the input array casted to unsigned.  If the
   array could not be created then returns NULL and errno is left as set by
   malloc().
**/
unsigned* a_ctou(char* S, unsigned n) {

  unsigned* r = malloc(n*sizeof(unsigned));
  if (!r)
    return 0;

  int i;
  for (i=0; i<n; i++)
    r[i] = (unsigned) S[i];

  return r;
}




/**
   \brief Load an array from a text file.

   a_load() loads n numbers of a given type from a text file into a new array
   and stores its address in *A. The memory region previously pointed to by *A,
   if any, will not be touched.

   The file format should be a blank separated succession of numbers of the
   given type.

   \param type The type of *A.  Valid types are 'i' (int) or 'd' (double).

   \return On success it returns 1 and sets *A to the address of the new array.
   On failure, it returns 0 and errno remains set as by fopen() on failure, or
   remains set as set by malloc() on failure or is set to EILSEQ if the input
   file has an incorrect format.
**/
int a_load(void** A, char type, unsigned n, char* filename) {

  FILE* f = fopen(filename,"rt");
  if (!f) return 0;

  void* V;

  if (type == 'i')
    V = malloc(sizeof(int)*n);
  else if (type == 'd')
    V = malloc(sizeof(double)*n);
  else
    V = 0;

  if (!V) {
    int err = errno;
    fclose(f);
    errno = err;
    return 0;
  }

  int i, k;
  for (i=0; i<n; i++) {
    if (type == 'i')
      k = fscanf(f,"%d ",((int*)V)+i);
    else if (type == 'd')
      k = fscanf(f,"%lf ",((double*)V)+i);
    else
      k = 0;

    if (k != 1) {
      free(V);
      fclose(f);
      errno = EILSEQ;
      return 0;
    }
  }

  fclose(f);
  *A = V;
  return 1;
}




/**
   \brief Print an array.

   a_print() prints each of the n trailing elements of A in a line followed by a
   space.

   \param type The type of *A.  Valid types are 'i' (int), 'u' (unsigned), 'h'
   (short), 'x' (u32), 'd' (double), c (char printed as int) or 's' (string).
**/
void a_print(void* A, char type, unsigned n) {

  if (A) {
    int i;

    for (i=0; i<n; i++) {
      if (type == 'i')
        printf("%d ",*(((int*)A)+i));
      else if (type == 'u')
        printf("%u ",*(((unsigned*)A)+i));
      else if (type == 'h')
        printf("%hd ",*(((short*)A)+i));
      else if (type == 'x')
        printf("%"PRIu32" ",*(((u32*)A)+i));
      else if (type == 'd')
        printf("%lf ",*(((double*)A)+i));
      else if (type == 'c')
        printf("%hhd ",*((char*)A+i));
      else
        printf("%s ",*((char**)A+i));
    }
  }
  printf("\n");
}




/**
   \brief Print an array in tabular form.

   a_print_as_matrix() prints the n*m elements trailing elements of array A in
   tabular layout.  Each element in a line is followed by a space.

   \param type The type of *A.
   Implemented types are 'i' (int) and 'd' (double).
**/
void a_print_as_matrix(void* A, char type, unsigned n, unsigned m) {

  int i,j;
  for (i=0; i<n; i++) {
    for (j=0; j<m; j++) {
      if (type == 'i')
        printf("%d ",*(((int*)A)+i*m+j));
      else if (type == 'd')
        printf("%f ",*(((double*)A)+i*m+j));
    }
    printf("\n");
  }
}




/**
   \brief Allocate a ragged array.

   m_alloc() allocates a n x m ragged array of a given type.

   \param type The type of array elements.
   Implemented types are 'i' (int), 'l' (long), 'd' (double), 'y' (u64).

   \return On success it returns the address of a new ragged array.  If the
   array could not be created then it returns 0 and errno remains as set by
   malloc().
**/
void** m_alloc(char type, unsigned n, unsigned m) {

  void** A;

  switch (type) {
  case 'd': A = malloc(n*sizeof(double*)); break;
  case 'i': A = malloc(n*sizeof(int*)); break;
  case 'l': A = malloc(n*sizeof(long*)); break;
  case 'y': A = malloc(n*sizeof(u64*)); break;
  default: A = 0;
  }

  if (!A)
    return 0;

  int i;
  for (i=0; i<n; i++) {

    switch (type) {
    case 'd': A[i] = malloc(m*sizeof(double)); break;
    case 'i': A[i] = malloc(m*sizeof(int)); break;
    case 'l': A[i] = malloc(m*sizeof(long)); break;
    case 'y': A[i] = malloc(m*sizeof(u64)); break;
    }

    if (!A[i]) {
      int err = errno;
      while (--i)
        free(A[i]);
      free(A);
      errno = err;
      return 0;
    }
  }

  return A;
}




/**
   \brief Allocate a ragged array and set it to zero.

   m_allocz() allocates a n x m ragged array of a given type filled with zeros.

   \param type The type of array elements.
   Implemented types are 'i'(int), 'l' (long), 'd'(double), 'y'(u64).

   \return On success it returns the address of the new ragged array.  If the
   array could not be created it then returns NULL and errno remains as set by
   calloc().
**/
void** m_allocz(char type, unsigned n, unsigned m) {

  void** A;

  switch (type) {
  case 'd': A = calloc(n,sizeof(double*)); break;
  case 'i': A = calloc(n,sizeof(int*)); break;
  case 'l': A = calloc(n,sizeof(long*)); break;
  case 'y': A = calloc(n,sizeof(u64*)); break;
  default: printf("No such type."); exit(1);
  }

  if (!A)
    return 0;

  int i;
  for (i=0; i<n; i++) {

    switch (type) {
    case 'd': A[i] = calloc(m,sizeof(double)); break;
    case 'i': A[i] = calloc(m,sizeof(int)); break;
    case 'l': A[i] = calloc(m,sizeof(long)); break;
    case 'y': A[i] = calloc(m,sizeof(u64)); break;
    }

    if (!A[i]) {
      int err = errno;
      while (--i)
        free(A[i]);
      free(A);
      errno = err;
      return 0;
    }
  }

  return A;
}




/**
  \brief Free a ragged array.

  m_free() frees a ragged array with n rows.
**/
void m_free(void** M, unsigned n) {

  if (!M) return;

  int i;
  for (i=0; i<n; i++)
    free(M[i]);
  free(M);
}




/**
   \brief Write a ragged array to a file.

   m_write() writes the submatrix [0,n-1],[0,m-1] of A to a text file with
   newline as the row separator and space as row separator.

   \param type The type of *A.
   Implemented types are 'i' (int) 'd' (double) 'y' (u64).

   \return On success it returns 1.  On failure it returns 0 and errno remains
   as set by fopen().
**/
int m_write(void** A, char type, unsigned n, unsigned m, char* filename) {

  return m_writes(A,type,n,m," ","\n",filename);
}




/**
   \brief Write a ragged array to a file.

   m_writes() writes the submatrix [0,n-1],[0,m-1] of A to a text file with
   given row and column separators.

   \param type The type of *A.
   Implemented types are 'i' (int) 'd' (double) 'y' (u64).

   \return On success it returns 1.  On failure it returns 0 and errno remains
   as set by fopen().
**/
int m_writes(void** A, char type, unsigned n, unsigned m,
             char* column_sep, char* row_sep, char* filename) {

  FILE* f = fopen(filename,"w");
  if (!f) return 0;

  int i,j;
  for (i=0; i<n; i++) {
    for (j=0; j<m; j++) {
      switch (type) {
      case 'i': fprintf(f,"%d%s",((int**)A)[i][j],column_sep); break;
      case 'd': fprintf(f,"%f%s",((double**)A)[i][j],column_sep); break;
      case 'y': fprintf(f,"%"PRIu64"%s",((u64**)A)[i][j],column_sep); break;
      }
      fprintf(f,"%s",row_sep);
    }
  }

  fclose(f);
  return 1;
}




/**
   \brief Load a matrix from a text file.

   m_load() loads n*m numbers of a given type from a text file into a new n x m
   ragged array and stores its address in ***M.  The memory previously pointed
   to by M will not be touched.

   The file should be a blank separated succession of numbers of the given type.

   \param type The type of **M.
   Implemented types are 'i' (int), 'd' (double), 'y' (uint64_t).

   \return On success returns it 1 and sets *M to the address of the new array.
   On failure it returns 0 and errno remains set either by fopen() or by
   malloc() on failure, or is set to EILSEQ if the input file has an unexpected
   format.
**/
int m_load(void*** M, char type, size_t n, size_t m, char* filename) {

  FILE* f = fopen(filename,"rt");
  if (!f)
    return 0;

  void **A = m_alloc(type,n,m);
  if (!A)
    return 0;

  int i,j,k;
  for (i=0; i<n; i++) {
    for (j=0; j<m; j++) {

      if (type == 'd')
        k = fscanf(f,"%lf",&((double**) A)[i][j]);
      else if (type == 'i')
        k = fscanf(f,"%d",&((int**) A)[i][j]);
      else if (type == 'y')
        k = fscanf(f,"%"SCNu64,&((u64**) A)[i][j]);
      else
        k = 0;

      if (k != 1) {
        fclose(f);
        m_free(A,n);
        errno = EILSEQ;
        return 0;
      }
    }
  }

  fclose(f);
  *M = A;
  return 1;
}



/**
   \brief Loads a lower triangular matrix from a text file into a symmetric matrix.

   This function loads numbers of a given type from a text file into a new n x n
   symmetric array.  The address of the new array is stored in ***M.  The memory
   previously pointed to by ***M will not be touched.

   The file should be organized as:

   <pre>
   M[1,0]
   M[2,0] M[2,1]
   M[3,0] M[3,1] M[3,2]
   ...
   M[n-1,0] ... M[n-1,n-2]
   </pre>

   \param type The type of **M.
   Implemented types are 'i' (int), 'd' (double), 'y' (uint64_t).

   \return On success returns 1 and sets *M to the address of the new array.  On
   failure, returns 0 and errno remains set either by fopen() or by malloc() on
   failure, or is set to EILSEQ if the input file has an incorrect format.
**/
int m_load_ltm(void*** M, char type, size_t n, char* filename) {

  FILE* f = fopen(filename,"rt");
  if (!f)
    return 0;

  void **A = m_allocz(type,n,n);
  if (!A)
    return 0;

  int i,j,k;
  for (i=1; i<n; i++) {
    for (j=0; j<i; j++) {

      if (type == 'd') {
        k = fscanf(f,"%lf",&((double**) A)[i][j]);
      }
      else if (type == 'i') {
        k = fscanf(f,"%d",&((int**) A)[i][j]);
      }
      else if (type == 'y') {
        k = fscanf(f,"%"SCNu64,&((u64**) A)[i][j]);
        ((u64**) A)[j][i] = ((u64**) A)[i][j];
      }
      else
        k = 0;

      if (k != 1) {
        fclose(f);
        m_free(A,n);
        errno = EILSEQ;
        return 0;
      }
    }
  }

  fclose(f);
  *M = A;
  return 1;
}




/**
   \brief Print a submatrix.

   This function prints elements in the submatrix M[frow,trow-1],[fcol,tcol-1].

   \param type The type of *M.
   Implemented types are 'i' (int) 'd' (double) 'y' (uint64_t).
**/
void m_prints(void** M, char type, size_t frow, size_t trow, size_t fcol, size_t tcol) {

  int i,j;

  for (i=frow; i<trow; i++) {
    //printf("%d ",i);

    for (j=fcol; j<tcol; j++) {

      if (type == 'i')
        printf("%d ",((int**) M)[i][j]);
      else if (type == 'd')
        printf("%.2lf ",((double**) M)[i][j]);
      else if (type == 'y')
        printf("%"PRIu64" ",((u64**) M)[i][j]);
    }

    printf("\n");
  }
}



/**
   \brief Print a matrix.

   \param type The type of *M.
   Implemented types are 'i' (int) 'd' (double) 'y' (uint64_t).
**/
void m_print(void** M, char type, size_t n, size_t m) {

  m_prints(M,type,0,n,0,m);
}



/**
   \brief Allocate a lower triangular matrix.

   This function allocates an order n lower triangular matrix of a given type.

   \param type The type of matrix elements.
   Implemented types are 'i' (int) 'd' (double).

   \return On success returns the address of the new ltm.  If the matrix could
   not be created then returns 0 and errno remains set by malloc() on failure.
**/
void** ltm_alloc(char type, unsigned n) {

  void** A;

  if (type == 'i')
    A = malloc(n*sizeof(int*));
  else if (type == 'd')
    A = malloc(n*sizeof(double*));
  else
    A = 0;

  if (!A)
    return 0;

  A[0] = 0;
  int i;
  for (i=1; i<n; i++) {
    if (type == 'i')
      A[i] = malloc(i*sizeof(int));
    else if (type == 'd')
      A[i] = malloc(i*sizeof(double));

    if (!A[i]) {
      int err = errno;
      while (--i)
        free(A[i]);
      free(A);
      errno = err;
      return 0;
    }
  }

  return A;
}



/**
   \brief Creates a copy of a stricly lower triangular matrix.

   \param S The source ltm.
   \param type The type of S elements.
   Implemented types are 'i' (int) 'd' (double).
   \param n The order of S.

   \return On success returns the address of the new ltm.  If the matrix could
   not be created then returns 0 and errno remains set by malloc() on failure.
**/
void** ltm_dup(void** S, char type, unsigned n) {

  void** A;

  switch (type) {
  case 'i': A = malloc(n*sizeof(int*)); break;
  case 'd': A = malloc(n*sizeof(double*)); break;
  default:  A = 0;
  }

  if (!A) return 0;

  int i;

  A[0] = 0;

  for (i=1; i<n; i++) {

    switch (type) {
    case 'i': A[i] = malloc(i*sizeof(int)); break;
    case 'd': A[i] = (double*) malloc(i*sizeof(double)); break;
    }

    if (!A[i]) {
      int err = errno;
      while (--i)
        free(A[i]);
      free(A);
      errno = err;
      return 0;
    }

    // int j;
    // for (j=0; j<i; j++) {
    //   switch (type) {
    //   case 'i':  break;
    //   case 'd': ((double**)A)[i][j] = ((double**)S)[i][j]; break;
    //   }
    // }

    switch (type) {
    case 'i': memcpy(((int**)A)[i],((int**)S)[i],i*sizeof(int)); break;
    case 'd': memcpy(((double**)A)[i],((double**)S)[i],i*sizeof(double)); break;
    }
  }
  return A;
}



/**
  \brief Releases a lower triangular matrix.

  This function releases a lower triangular matrix with n rows.
**/
void ltm_free(void** M, unsigned n) {
  m_free(M,n);
}


/**
   \brief Load a lower triangular matrix from a text file.

   This function loads (n-1)*(n-2)/2 numbers from a text file into a new order n
   lower triangular matrix pointed and stores its address in *M.

   The memory region previously pointed to by M, if any, will not be touched.

   The file format should be a blank separated succession of numbers of the
   given type.

   \param type The type of **M.  Valid types are 'i' (int) or 'd' (double).

   \return On success returns 1 and stores the address of the new ltm in *M.
   On failure returns 0 and errno remains set either by fopen() on failure or by
   malloc() on failure or by fscanf() on failure.
**/
int ltm_load(void*** M, char type, unsigned n, char* filename) {

  FILE* f = fopen(filename,"rt");
  if (!f)
    return 0;

  void **A = ltm_alloc(type,n);
  if (!A)
    return 0;

  int i,j,k;
  for (i=1; i<n; i++) {
    for (j=0; j<i; j++) {

      if (type == 'd')
        k = fscanf(f,"%lf",&((double**) A)[i][j]);
      else if (type == 'i')
        k = fscanf(f,"%d",&(((int**) A))[i][j]);
      else
        k = 0;

      if (k != 1) {
        int err = errno;
        fclose(f);
        m_free(A,n);
        errno = err;
        return 0;
      }
    }
  }

  fclose(f);
  *M = A;
  return 1;
}




/**
   \brief Print a lower triangular matrix.

   This function prints an order n lower triangular matrix of a given type.

   \param type The type of *A.
   Implemented types are 'i' (int) 'd' (double) 'y' (u64).
**/
void ltm_print(void** A, char type, unsigned n) {

  int i,j;

  for (i=0; i<n-1; i++) {
    printf("%d ",i);
  }

  printf("\n");
  //printf("null\n");
  for (i=1; i<n; i++) {
    printf("%d ",i);
    for (j=0; j<i; j++) {
      switch (type) {
      case 'i': printf("%d ",((int**)A)[i][j]); break;
      case 'd': printf("%.2f ",((double**)A)[i][j]); break;
      case 'y': printf("%"PRIu64" ",((u64**)A)[i][j]); break;
      }
    }
    printf("\n");
  }
}



/**
   \brief Write a lower triangular matrix to file.

   This function writes an order n lower triangular matrix of a given type to a
   file.

   \param type The type of *A.
   Implemented types are 'i' (int) 'd' (double) 'y' (u64).

   \return On success it returns 1.  On failure it returns 0 and errno remains
   as set by fopen().
**/
int ltm_write(void** A, char type, unsigned n, char* filename) {

  FILE* f = fopen(filename,"w");
  if (!f) return 0;

  int i,j;
  //printf("null\n");
  for (i=1; i<n; i++) {
    for (j=0; j<i; j++) {
      switch (type) {
      case 'i': fprintf(f,"%d ",((int**)A)[i][j]); break;
      case 'd': fprintf(f,"%.2f ",((double**)A)[i][j]); break;
      case 'y': fprintf(f,"%"PRIu64" ",((u64**)A)[i][j]); break;
      }
    }
    fprintf(f,"\n");
  }

  fclose(f);
  return 1;
}



/**
   \brief Creates a binary string from an unsigned integral.
*/
char* u8toa(u8 a, char* buffer) {

  buffer += 8;
  *buffer-- = 0;

  int i;
  for (i = 7; i >= 0; i--) {
    *buffer-- = (a & 1) + '0';
    a >>= 1;
  }

  return buffer+1;
}




/**
@copydoc u8toa
*/
char* u64toa(u64 a, char* buffer) {

  buffer += 64;
  *buffer-- = 0;

  int i;
  for (i = 63; i >= 0; i--) {
    *buffer-- = (a & 1) + '0';
    a >>= 1;
  }

  return buffer+1;
}


/**
@copydoc u8toa
*/
char* u32toa(u32 a, char* buffer) {

  buffer += 32;
  *buffer-- = 0;

  int i;
  for (i = 31; i >= 0; i--) {
    *buffer-- = (a & 1) + '0';
    a >>= 1;
  }

  return buffer+1;
}






// /**
// Reads 8-byte doubles from a file into a lower triangular matrix (as a
//    ragged array).
//
//    Space for M is allocated by the function, and the memory region previously
//    pointed to by M, if any, will not be touched.
//
//    The matrix dimension d will will be calculated from file length, that should
//    be equal to d(d-1)/2.
//
//    This function returns the matrix dimension.  On failure, it returns -1 and
//    sets errno as follows.
//
//    - The same set by stat() or fopen().
//    - ENOMEM if not enough memory is available.
//    - EILSEQ if the function fails to read the matrix completely.
//    - EOVERFLOW if the size of a double is not 8 bytes.
// **/
//
// int ltmd_load(char* filename, double ***M) {
//
//   if (sizeof(double) != 8) {
//     errno = EOVERFLOW;
//     return -1;
//   }
//
//   struct stat s;
//   if (stat(filename, &s) < 0)
//     return 0;
//
//   n = s.st_size / sizeof(double);
//   d = (1+sqrt(1+sizeof(double)*n))/2;
//
//   FILE* f = fopen(filename,"rb");
//   if (!f)
//     return -1;
//
//   double** A = ltmd_alloc(d);
//   if (!A) {
//     errno = ENOMEM;
//     return -1;
//   }
//
//   int i;
//
//   A[0] = 0;
//   for (i=1; i<d; i++) {
//     if (fread(A[i],sizeof(double),i,f) != i) {
//       m_free((void**)A,d);
//       errno = EILSEQ;
//       return -1;
//     }
//   }
//
//   fclose(f);
//
//   *M = A;
//   return d;
// }
// /**
//    Writes a sequence of d(d-1)/2 random doubles in [0,1) to a file.  The sequence
//    corresponds to a lower triangular matrix of order d.
//
//    Returns d on success.  On failure, returns -1 and sets errno to indicate the
//    error, EOVERFLOW meaning that the double size is not 8 bytes.
// **/
//
// int ltmd_write_random(int d, char* matf, unsigned int seed) {
//
//   if (sizeof(double) != 8) {
//     errno = EOVERFLOW;
//     return -1;
//   }
//
//   long long n = ((long long) d)*(d-1)/2;
//
//   FILE* f = fopen(matf,"wb");
//   if (!f)
//     return -1;
//
//
//   long long i;
//   for (i=0; i<n; i++) {
//     double r = (double) random()/(RAND_MAX+1.0);
//     if (fwrite(&r,sizeof(double),1,f) != 1)
//       return -1;
//   }
//
//   fclose(f);
//
//   return d;
// }
