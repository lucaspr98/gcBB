/*
  Guilherme P. Telles.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#include "distance.h"
#include "abbrevs.h"
#include "arrays.h"
#include "estringue.h"
#include "utils.h"


/**
   \brief Cosine distance.

   Evaluate cosine distances among n vectors with d dimensions stored linearly
   in F and store in D.

   The first vector in F is F[0],...,F[d-1], the second is F[d],...,F[2d-1] and
   so on.

   If a non-null index I on F is provided, then the first vector will be
   F[I[0]d],...,F[I[0](d+1)-1], the second will be F[I[1]d],...,F[I[1](d+1)-1]
   and so on.

   \param F n vectors with d dimensions.
   \param n The number of vectors.
   \param d The number of dimensions of each vector.
   \param I An index on F or NULL.
   \param D An existing strictly lower triangular matrix of order n (see arrays.h).

   \return On success it returns 1.  On failure it returns 0 and either errno
   remains as set by malloc() or errno is set to EOVERFLOW if any distance
   didn't fit in a double.
**/
int cosine_distance(double* F, size_t n, size_t d, int* I, double** D) {

  double* en = malloc(n*sizeof(double));
  if (!en) return 0;

  int i,j,t;

  for (i=0; i<n; i++) {
    double* Fi = F+(I?I[i]:i)*d;
    double aux = 0;
    for (t=0; t<d; t++)
      aux += *(Fi+t) * *(Fi+t);
    en[i] = sqrt(aux);
  }

  for (i=0; i<n; i++) {
    double* Fi = F+(I?I[i]:i)*d;
    for (j=0; j<i; j++) {
      double* Fj = F+(I?I[j]:j)*d;
      double dij = 0;
      for (t=0; t<d; t++)
        dij += *(Fi+t) * *(Fj+t);

      double np = en[i]*en[j];
      D[i][j] = 1 - (np == 0 ? -1 : dij/np);

      if (isinf(D[i][j])) {
        free(en);
        errno = EOVERFLOW;
        return 0;
      }
    }
  }

  free(en);
  return 1;
}



/**
   \brief Euclidean distance.

   Evaluate Euclidean distances among n vectors with d dimensions stored linearly
   in F and store in D.

   The first vector in F is F[0],...,F[d-1], the second is F[d],...,F[2d-1] and
   so on.

   If a non-null index I is provided, then the first vector will be
   F[I[0]d],...,F[I[0](d+1)-1], the second will be F[I[1]d],...,F[I[1](d+1)-1]
   and so on.

   \param F n vectors with d dimensions.
   \param n The number of vectors.
   \param d The number of dimensions of each vector.
   \param I An index on F or NULL.
   \param D An existing strictly lower triangular matrix of order n (see arrays.h).

   \return On success it returns 1.  On failure it returns 0 and errno is set to
   EOVERFLOW if any distance didn't fit in a double.
**/
int euclidean_distance(double* F, size_t n, size_t d, int* I, double** D) {

  int i,j,t;

  for (i=0; i<n; i++) {
    double* Fi = F+(I?I[i]:i)*d;
    for (j=0; j<i; j++) {
      double* Fj = F+(I?I[j]:j)*d;
      double dij = 0;
      for (t=0; t<d; t++) {
        double aux = *(Fi+t) - *(Fj+t);
        dij += aux * aux;
      }

      D[i][j] = sqrt(dij);

      if (isnan(D[i][j])) {
        errno = EOVERFLOW;
        return 0;
      }
    }
  }

  return 1;
}



/**
   \brief Allocate a new dmat structure.

   Create a new dmat structure and initialize its fields as follows: n is set
   to the homonym argument value, M is set to a new, non-initialized, order n
   lower triangular matrix, and labels is set to a new array of n null pointers.

   \return On success it returns the address of a new dmat structure.  On failure
   returns NULL and errno remains set as bu malloc() on failure.
**/
dmat* dmat_alloc(int n) {

  int err;

  dmat* R = malloc(sizeof(dmat));
  if (!R) goto ENOMEMH;

  R->n = n;

  R->labels = (char**) calloc(n,sizeof(char*));
  if (!R->labels) goto ENOMEMH;

  R->M = (double**) ltm_alloc('d',n);
  if (!R->M) goto ENOMEMH;

  return R;

 ENOMEMH:
  err = errno;
  dmat_free(R);
  errno = err;

  return 0;
}



/**
   \brief Release a dmat structure and its fields.
**/
void dmat_free(dmat* T) {

  if (!T)
    return;
  if (T->labels)
    m_free((void**)T->labels,T->n);
  if (T->M)
    ltm_free((void**)T->M,T->n);
  free(T);
}



/**
   \brief Copy a dmat structure.

   Create a through copy of a dmat structure.

   \return On success it returns the address of a new dmat structure.  On failure
   returns NULL and errno remains set as bu malloc() on failure.
**/
dmat* dmat_copy(dmat* D) {
  int i;
  int err;

  dmat* R = malloc(sizeof(dmat));
  if (!R) goto ENOMEMH;

  R->n = D->n;

  if (D->labels) {
    R->labels = (char**) calloc(D->n,sizeof(char*));
    if (!R->labels) goto ENOMEMH;

    for (i=0; i<D->n; i++)
      if (D->labels[i]) {
        R->labels[i] = strdup(D->labels[i]);
        if (!R->labels[i]) goto ENOMEMH;
      }
  }
  else
    R->labels = NULL;

  R->M = (double**) ltm_alloc('d',D->n);
  if (!R->M) goto ENOMEMH;

  for (i=1; i<D->n; i++)
    memcpy(R->M[i],D->M[i],i*sizeof(double));

  return R;

 ENOMEMH:
  err = errno;
  dmat_free(R);
  errno = err;

  return 0;
}



/**
  \brief Read a dmat file.

  This function loads distances from a file formatted as the example below into
  a dmat structure.

  The labels section is optional.  If present, there must be at least n words
  separated by blanks in the section; the first n words will be read as labels.
  If absent, the field labels in the returning dmat will be NULL.

  There may be excess numbers in the distances section.  The first (size^2 -
  size)/2 floating-point numbers will be read.

  <pre>
  [size]
  4
  [labels]
  a b c d
  [distances]
  D(1,0)
  D(2,0) D(2,1)
  D(3,0) D(3,1) D(3,2)
  </pre>

  \return On success it returns a new dmat structure.  On failure it returns
  NULL and either errno remains set as by fopen() or malloc() on failure or
  errno is set to EILSEQ to indicated that parsing the file failed.  File format
  error checking is not comprehensive, the o
**/
dmat* dmat_read(char* filename) {

  FILE* f = fopen(filename,"r");
  if (!f) return 0;

  int i,j,err;
  dmat* R = 0;

  char* word = malloc(100*sizeof(char));
  if (!word) goto ENOMEMH;

  int n;
  if (fscanf(f," [size] %d ",&n) != 1) goto EILSEQH;

  R = dmat_alloc(n);
  if (!R) goto ENOMEMH;

  if (fscanf(f,"%s ",word) != 1) goto EILSEQH;

  if (!strcmp(word,"[labels]")) {
    for (i=0; i<n; i++) {
      if (fscanf(f,"%s ",word) != 1) goto EILSEQH;
      R->labels[i] = strdup(word);
      if (!R->labels[i]) goto ENOMEMH;
    }
    if (fscanf(f,"%s ",word) != 1) goto EILSEQH;
  }
  else {
    free(R->labels);
    R->labels = NULL;
  }

  while (strcmp(word,"[distances]"))
    if (fscanf(f,"%s ",word) != 1) goto EILSEQH;

  for (i=1; i<n; i++)
    for (j=0; j<i; j++)
      if (fscanf(f,"%lf ",&(R->M[i][j])) != 1) goto EILSEQH;

  fclose(f);
  free(word);
  return R;

 ENOMEMH: err = errno; goto STDH;
 EILSEQH: err = EILSEQ; goto STDH;
 STDH:
  fclose(f);
  dmat_free(R);
  free(word);
  errno = err;
  return 0;
}



/**
  \brief Write a dmat to file.

  See dmat_read() for a file format description.

  \param D A dmat.
  \param cast_to_long If nonzero distances are cast to long and written to file as integral values.
  \param filename The output file name.

  \return On success it returns 1.  On failure it returns 0 and errno remains
  as set by fopen().
**/
int dmat_write(dmat* D, int cast_to_long, char* filename) {
  int i, j;

  FILE* f = fopen(filename,"w");
  if (!f) return 0;

  fprintf(f,"[size]\n%d\n",D->n);

  if (D->labels) {
    fprintf(f,"[labels]\n");
    for (i=0; i<D->n; i++) {
      fprintf(f,"%s ",D->labels[i]);
    }
  }

  fprintf(f,"[distances]\n");

  for (i=1; i<D->n; i++) {
    for (j=0; j<i; j++) {
      if (cast_to_long)
        fprintf(f,"%ld ",(long) D->M[i][j]);
      else
        fprintf(f,"%.16lf ",D->M[i][j]);
    }
    fprintf(f,"\n");
  }

  fclose(f);
  return 1;
}



/**
  \brief Write symmetric integral distances to a dmat file.

  \param M A symmetric matrix with size n times n, at least.
  \param n The order of M.
  \param filename The file name.

  See dmat_read() for a file format description.

  \return On success it returns 1.  On failure it returns 0 and errno remains
  as set by fopen().
**/
int write_as_dmat(int** M, int n, char* filename) {
  int i,j;
  FILE* f = fopen(filename,"w");
  if (!f) return 0;

  fprintf(f,"[size]\n%d\n",n);

  fprintf(f,"[distances]\n");

  for (i=1; i<n; i++) {
    for(j=0; j<i; j++) {
      fprintf(f,"%d ",M[i][j]);
    }
    fprintf(f,"\n");
  }

  fclose(f);
  return 1;
}
