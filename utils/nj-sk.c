/**
  Guilherme P. Telles.
**/

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <errno.h>
#include <math.h>
#include <string.h>

#include "arrays.h"
#include "graph.h"



/**
   \brief Neighbor-Joining with Studier and Keppler's equations.

   Neighbor-Joining implementing the equations proposed by Studier and Keppler
   (J.A. Studier and K.J. Keppler. A note on the Neighbor-Joining algorithm of
   Saitou and Nei.  Mol. Biol. Evol. v.5, 1988.)

   \param D An order n strictly lower triangular matrix with distances among
   OTUs.  Its contents will not be preserved, but it won't be reallocated or
   freed.

   \param n The number of OTUs.

   \returns An undirected acyclic graph where OTUs are the vertices with labels
   in [0,n-1] and hypothetical ancestors are the vertices with labels in
   [n,2n-3].  In the graph OTUs have their flag field set to 1 and ancestors
   have their flag field set to 0.  On failure this function returns NULL and
   sets errno to ENOMEM.
**/
graph* nj_sk(double** D, unsigned n) {

  int i,j,k;
  graph* T = NULL;
  int* map = NULL;
  double* R = NULL;

  if (n == 1) {
    T = gr_alloc('u',1);
    if (!T) goto NOMEMH;
    T->V[0].flag = 1;
    return T;
  }

  if (n == 2) {
    T = gr_alloc('u',2);
    if (!T) goto NOMEMH;
    T->V[0].flag = 1;
    T->V[1].flag = 1;
    if (!gr_add_edge(T,0,1,D[1][0]))
      goto NOMEMH;
    return T;
  }


  // The tree:
  T = gr_alloc('u',2*n-2);
  if (!T) goto NOMEMH;

  // Store the index of the OTU at row i of D in map[i].  Initially map[i] = i.
  // Set flag in OTU vertices to 1.
  map = malloc(n*sizeof(int));
  if (!map) goto NOMEMH;

  for (i=0; i<n; i++) {
    map[i] = i;
    T->V[i].flag = 1;
  }

  // Evaluate the sum of distances from one OTU to the others:
  R = (double*) malloc(n*sizeof(double));
  if (!R) goto NOMEMH;

  for (k=0; k<n; k++) {
    R[k] = 0;
    for (i=0; i<k; i++)
      R[k] += D[k][i];
    for (i=k+1; i<n; i++)
      R[k] += D[i][k];
  }

  // The next hypothetical node:
  int hyp = n;

  while (n > 3) {

    // Evaluate Q[i,j] and get minimum:
    double q, qmin = DBL_MAX;
    int imin = 0, jmin = 0;

    for (i=1; i<n; i++) {
      for (j=0; j<i; j++) {
        q =  ((n-2) * D[i][j]) - R[i] - R[j];
        if (q < qmin) {
          qmin = q;
          imin = i;
          jmin = j;
        }
      }
    }

    // Branches ik and jk:
    double lik = ((n-2)*D[imin][jmin] + R[imin] - R[jmin]) / ((double)(2*n-4));
    double ljk = D[imin][jmin] - lik;

    if (!gr_add_edge(T,map[imin],hyp,lik) || !gr_add_edge(T,map[jmin],hyp,ljk))
      goto NOMEMH;

    // Remove jmin and imin from R:
    for (k=0; k<jmin; k++)
      R[k] -= D[jmin][k];
    for (k=jmin+1; k<n; k++)
      R[k] -= D[k][jmin];

    for (k=0; k<imin; k++)
      R[k] -= D[imin][k];
    for (k=imin+1; k<n; k++)
      R[k] -= D[k][imin];

    // The new node goes on jmin:
    map[jmin] = hyp;

    R[jmin] = 0;
    for (k=0; k<jmin; k++) {
      D[jmin][k] = (D[jmin][k] + (k<imin ? D[imin][k] : D[k][imin]) - D[imin][jmin]) / 2;
      R[jmin] += D[jmin][k];
      R[k] += D[jmin][k];
    }

    for (k=jmin+1; k<n; k++) {
      if (k != imin) {
        D[k][jmin] = (D[k][jmin] + (k<imin ? D[imin][k] : D[k][imin]) - D[imin][jmin]) / 2;
        R[jmin] += D[k][jmin];
        R[k] += D[k][jmin];
      }
    }

    // Move n-1 on imin:
    if (n-1 != imin) {
      for (k=0; k<imin; k++)
        D[imin][k] = D[n-1][k];
      for (k=imin+1; k<n-1; k++)
        D[k][imin] = D[n-1][k];

      R[imin] = R[n-1];
      map[imin] = map[n-1];
    }

    hyp++;
    n--;
  }

  // 3 points:
  double x = (D[1][0]+D[2][0]-D[2][1])/2;
  double y = (D[1][0]+D[2][1]-D[2][0])/2;
  double z = (D[2][0]+D[2][1]-D[1][0])/2;

  if (!gr_add_edge(T,map[0],hyp,x) ||
      !gr_add_edge(T,map[1],hyp,y) ||
      !gr_add_edge(T,map[2],hyp,z))
    goto NOMEMH;

  free(map);
  free(R);

  return T;

 NOMEMH:
  gr_free(T);
  free(map);
  free(R);
  errno = ENOMEM;
  return 0;
}
