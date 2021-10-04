/**
   \file distance.h
   \brief Functions for distances.

   \internal Guilherme P. Telles, June 2017, sep 2019.
**/

#ifndef DISTANCEH
#define DISTANCEH

/**
   \brief A symmetric distance matrix.  The matrix is a lower triangular matrix
   of order n represented as a ragged array with the following form:

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
**/
struct dmat {
  int n;     ///<\brief The matrix order.
  double** M;     ///<\brief A lower triangular matrix with order n.
  char** labels;  ///<\brief An array of n labels.
};

typedef struct dmat dmat;


int cosine_distance(double* F, size_t n, size_t d, int* I, double** D);
int euclidean_distance(double* F, size_t n, size_t d, int* I, double** D);

dmat* dmat_alloc(int n);
void dmat_free(dmat* T);

dmat* dmat_copy(dmat* D);

dmat* dmat_read(char* filename);
int dmat_write(dmat* D, int cat_to_long, char* filename);

int write_as_dmat(int** M, int n, char* filename);

#endif
