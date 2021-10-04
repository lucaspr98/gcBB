// Guilherme P. Telles, 2020.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <getopt.h>

#include "arrays.h"
#include "distance.h"
#include "graph.h"
#include "utils.h"

#include "nj-sk.h"

#if USEMC==1
#include "malloc_count/malloc_count.h"
#endif


void help() {
  printf("Usage: nj-main -i input-dmat -d output-dot -n output-newich -m output-dmat \n");
  exit(1);
}


int main (int argc, char **argv) {

  char* in_dmat = NULL;
  char* out_dot = NULL;
  char* out_nw = NULL;
  char* out_dmat = NULL;

  extern int opterr;
  opterr = 0;

  extern int optind;
  int c;
  while ((c = getopt(argc, argv, "i:d:n:m:")) != -1) {

    switch (c) {

    case 'i':
      in_dmat = strdup(optarg);
      break;

    case 'd':
      out_dot = strdup(optarg);
      break;

    case 'n':
      out_nw = strdup(optarg);
      break;

    case 'm':
      out_dmat = strdup(optarg);
      break;

    default:
      help();
    }
  }

  if (optind < argc || optind == 1)
    help();

  dmat* D = dmat_read(in_dmat);
  if (!D) die("Unable to load %s.\n",in_dmat);

  clock_t t = 0;
  graph* T = NULL;

  t = clock();
  T = nj_sk(D->M,D->n);
  t = clock() - t;

  printf("time %Lf\n",(long double)t/CLOCKS_PER_SEC);

  if (out_dot) {
    if (!gr_write_dot(T,D->labels,D->n,1,in_dmat,out_dot))
      die("Unable to write %s.",out_dot);
  }

  if (out_nw) {
    if (!tree_write_newicks(T,D->labels,D->n,out_nw))
      die("Unable to write %s.",out_nw);
  }

  if (out_dmat) {
    ltm_free((void**)D->M,D->n);
    D->M = tree_apw(T,D->n);
    free(D->labels);
    D->labels = NULL;

    if (!dmat_write(D,0,out_dmat))
      die("Unable to write %s.",out_dmat);
  }

  dmat_free(D);
  gr_free(T);

  return 0;
}
