/**
   \file graph.h
   \brief A weighted graph, either oriented or non-oriented.

   A simple graph as adjacency lists.  Vertices are indexed from 0 to n-1.

   \internal Guilherme P. Telles.
**/


#ifndef GRAPHH
#define GRAPHH

//#include "abbrevs.h"
#include "list.h"


/**
   \brief An edge.
**/
struct edge {
  unsigned term; ///<\brief Terminal vertex index.
  double w;      ///<\brief Weight.
  int flag;      ///<\brief A general purpose flag.
  struct edge *next;
};

typedef struct edge edge;


/**
   \brief A vertex.
**/
struct vertex {
  edge* N;   ///<\brief The vertex's neighborhood.
  int flag;  ///<\brief A general purpose flag.
};

typedef struct vertex vertex;


/**
   \brief A graph.
**/
struct graph {
  int n;      ///<\brief The number of vertices.
  int m;      ///<\brief The number of edges.
  vertex *V;  ///<\brief The vertices, |V|=n.
  char type;  ///<\brief The graph type: 'd' for a directed graph, 'u' for an undirected graph.
  int flag;   ///<\brief A general purpose flag.
};

typedef struct graph graph;


graph* gr_alloc(char type, int n);
void gr_free(graph *G);

graph* gr_copy(graph* G);

int gr_add_edge(graph *G, int u, int v, double w);
int gr_remove_edge(graph *G, int u, int v);
int gr_remove_edges(graph *G, int u);

int gr_internal_edges(graph *G, int both_directions, int** edges, int* m);

void gr_neighbors(graph *G, int u, int* N, int* k);


void gr_exchange_indexes(graph* G, int u, int v);

edge* gr_edge(graph *G, int u, int v);
double gr_edge_weight(graph *G, int u, int v);

graph* gr_read(char* filename);
graph* gr_read_dimacs_uu(char* filename);

int gr_print(graph* G);
int gr_write_dot(graph* G, char** labels, int l, int print_edge_weights, char* title, char* file);
int gr_write_dot_images(graph* G, char** image_files, int l, int print_edge_weights, char* file);
int gr_write_json(graph* G, int print_edge_weights, char* file);

void gr_max_degree(graph* G, int* maxdeg, int* v);
void gr_eval_degrees(graph* G, int* maxdeg, int* v);

double** gr_apspw(graph* G);
int** gr_apd(graph* G, int k);

int gr_bfs(graph* G, int s);
int gr_diameter(graph* G);

list* gr_max_cliques1(graph* G);
list* gr_max_cliques2(graph* G);


char* tree_newick(graph *T);
char* tree_newicks(graph *T, char** labels, int l);
int tree_write_newicks(graph *T, char** labels, int l, char* file);

int tree_sspw(graph *G, int s, int* pi, double* weights);
double** tree_apw(graph* G, int k);

int tree_center(graph *G, int* c1, int* c2, int *dia);
int tree_nearest_leaf(graph *T, int s, int* l, int* d);

#endif
