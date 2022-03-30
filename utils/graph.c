/*
  Guilherme P. Telles, ?-2020.
*/

#define _GNU_SOURCE

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include "utils.h"
#include "graph.h"
#include "arrays.h"
#include "dynarray.h"
#include "list.h"



/**
  \brief Allocate a graph.

  \param type The graph type, 'd' for directed graph or 'u' for undirected graph.
  \param n The number of vertices.

  \return On success gr_alloc() returns the address of a new graph.  Otherwise
  it returns NULL and errno is left as set by malloc() on failure.
**/
graph* gr_alloc(char type, int n) {

  graph *G;
  if (!(G = malloc(sizeof(graph))))
    return 0;

  G->V = calloc(n,sizeof(vertex));
  if (!G->V && n>0) {
    int err = errno;
    free(G);
    errno = err;
    return 0;
  }

  G->type = type;
  G->n = n;
  G->m = 0;

  return G;
}



/**
   \brief Release a graph.
**/
void gr_free(graph* G) {

  if (!G) return;

  int i;
  for (i=0; i<G->n; i++) {
    edge* e = G->V[i].N;
    while (e) {
      edge* f = e;
      e = e->next;
      free(f);
    }
  }

  free(G->V);
  free(G);
}

/**
   Add adjacency uv to the beginning of adjacency list G[u].  It doesn't update G->m.
   Return 0 on malloc failure, 1 otherwise.
*/
int add_uv(graph* G, int u, int v, double w) {

  edge *uv = malloc(sizeof(edge));
  if (!uv) return 0;

  uv->term = v;
  uv->w = w;
  uv->flag = 0;
  uv->next = G->V[u].N;
  G->V[u].N = uv;
  return 1;
}




graph* gr_copy(graph* G) {
  int u;
  graph* H = gr_alloc(G->type,G->n);
  if (!H) return NULL;

  for (u=0; u<G->n; u++) {
    edge* e = calloc(1, sizeof(edge));
    for (e = G->V[u].N; e; e = e->next) {
      int st = add_uv(H,u,e->term,e->w);
      if (!st) {
        gr_free(H);
        return NULL;
      }
    }
  }

  H->m = G->m;

  for (u=0; u<G->n; u++) {
    H->V[u].flag = G->V[u].flag;
    edge *e, *f;
    for (e = G->V[u].N, f = H->V[u].N; e; e = e->next, f = f->next) {
      e->flag = f->flag;
    }
  }

  return H;
}






/**
   Remove adjacency uv from G.  It doesn't update G->m.
   Return 1 if the adjacency was removed, 0 if it didn't exist.
**/
int remove_uv(graph* G, int u, int v) {

  edge* p = 0;
  edge* e = G->V[u].N;

  while (e) {
    if (e->term == v) {
      if (p)
        p->next = e->next;
      else
        G->V[u].N = e->next;

      free(e);
      return 1;
    }

    p = e;
    e = e->next;
  }

  return 0;
}




/**
  \brief Add a weighted edge to a graph.

  If the graph is undirected then this function adds both uv and vu adjacencies.
  If u>=n or v>=n then the graph will be resized to have max(u,v)+1 vertices.

  \return If successful it returns 1. If an error occurs it returns 0 and errno
  remains as set by malloc() or realloc() on failure.
**/
int gr_add_edge(graph *G, int u, int v, double w) {

  // Resize V if necessary:
  if (u >= G->n || v >= G->n) {
    int nn = max(u,v)+1;
    vertex* W = realloc(G->V,nn*sizeof(vertex));
    if (!W)
      return 0;

    memset(W+G->n,0,(nn-G->n)*sizeof(vertex));
    G->V = W;
    G->n = nn;
  }

  // Add edge:
  int st = add_uv(G,u,v,w);
  if (!st)
    return 0;

  if (G->type == 'u') {
    st = add_uv(G,v,u,w);
    if (!st) {
      st = errno;
      remove_uv(G,u,v);
      errno = st;
      return 0;
    }
  }

  G->m++;
  return 1;
}



/**
  \brief Remove an edge from a graph.

  If the graph is undirected then this function removes both adjacencies uv and vu.

  \return The number of removed edges.
**/
int gr_remove_edge(graph *G, int u, int v) {

  int r = remove_uv(G,u,v);

  if (r) {
    if (G->type == 'u')
      remove_uv(G,v,u);
    G->m--;
  }

  return r;
}




/**
  \brief Remove every edge incident to a vertex.

  \return The number of removed edges.
**/
int gr_remove_edges(graph *G, int u) {

  int r = 0;

  while (G->V[u].N) {
    if (G->type == 'u')
      remove_uv(G,G->V[u].N->term,u);
    edge* e = G->V[u].N;
    G->V[u].N = G->V[u].N->next;
    free(e);
    G->m--;
    r++;
  }

  return r;
}


/**
   \brief Get the neighbors (out-neighbors) of v in G.

   This function does not modify graph flags.

   \param G The graph.
   \param v The vertex.
   \param maxdeg An existing array where vertex indices will be written.
   \param k The number of neighbors.
**/
void gr_neighbors(graph *G, int u, int* N, int* k) {

  *k = 0;
  edge* e = G->V[u].N;
  while (e) {
    N[*k] = e->term;
    (*k)++;
    e = e->next;
  }
}







/**
   \brief Swap the indexes of vertices u and v.

   Modify the graph adjusting neigborhoods and flags of u and v such that their
   indexes are exchanged.

   Time: O(D2(u)+D2(v)), where D2(.) is the number of vertices at distance 2, for
   undirected graphs. O(n+m) for directed graphs.
   Space: O(1).
**/
void gr_exchange_indexes(graph* G, int u, int v) {

  if (u == v) return;

  edge *e, *f;
  int i;

  if (G->type == 'u') {
    for (e = G->V[u].N; e; e = e->next)
      if (e->term == v)
        e->term = u;
      else {
        for (f = G->V[e->term].N; f->term != u; f = f->next) ;
        f->term = v;
      }

    for (e = G->V[v].N; e; e = e->next)
      if (e->term == u)
        e->term = v;
      else {
        for (f = G->V[e->term].N; f->term != v; f = f->next) ;
        f->term = u;
      }
  }
  else {
    for (i=0; i<G->n; i++)
      for (e = G->V[i].N; e; e = e->next)
        if (e->term == u)
          e->term = v;
        else if (e->term == v)
          e->term = u;
  }

  e = G->V[u].N;
  i = G->V[u].flag;

  G->V[u].N = G->V[v].N;
  G->V[u].flag = G->V[v].flag;

  G->V[v].N = e;
  G->V[v].flag = i;
}



/**
  \brief Retrieve a pointer to an edge.

  \return A pointer to uv or NULL.
**/
edge* gr_edge(graph *G, int u, int v) {

  edge* e = G->V[u].N;
  while (e) {
    if (e->term == v)
      break;
    e = e->next;
  }

  return e;
}



/**
  \brief Retrieve the weight of an edge.

  \return The weight of uv or DBL_MIN.
**/
double gr_edge_weight(graph *G, int u, int v) {

  edge* e = G->V[u].N;
  while (e) {
    if (e->term == v)
      break;
    e = e->next;
  }

  return e ? e->w : DBL_MIN;
}



/**
   \brief Read a graph from a file.

   Read a graph from a file with the first three lines with the format:

   <pre>
   vertices: n
   type: u/d
   edge weights: y/n
   </pre>

   The following lines have edges.  For a graph without edge weigths, each line
   must have two int integers. For a graph with edge weigths, each
   line must have two int integers and a double.

  \return On success it returns a new graph.  On failure it returns NULL and
  either errno remains set as by fopen() or malloc() on failure or errno is set
  to EBFONT to indicated that parsing the file failed.

**/
graph* gr_read(char *filename) {

  FILE *in;
  graph *G = NULL;
  int n, u, v, st;
  char type, weighted;

  if (!(in = fopen(filename,"r")))
    return 0;

  st = fscanf(in,"vertices: %d ",&n);
  if (st != 1) goto fmth;

  st = fscanf(in,"type: %c ",&type);
  if (st != 1) goto fmth;
  if (type != 'u' && type != 'd') goto fmth;

  st = fscanf(in,"edge weights: %c ",&weighted);
  if (st != 1) goto fmth;
  if (weighted != 'y' && weighted != 'n') goto fmth;

  G = gr_alloc(type,n);
  if (!G) goto nomemh;

  if (weighted == 'n') {
    while ((st = fscanf(in,"%d %d",&u,&v)) != EOF) {
      if (st != 2) goto fmth;
      gr_add_edge(G,u,v,0.0);
    }
  }
  else {
    double w;
    while ((st = fscanf(in,"%d %d %lf",&u,&v,&w)) != EOF) {
      if (st != 3) goto fmth;
      gr_add_edge(G,u,v,w);
    }
  }

  if (ferror(in)) {
    int err = errno;
    gr_free(G);
    errno = err;
    return NULL;
  }

  fclose(in);
  return G;


 fmth:
  errno = EBFONT;

 nomemh:
  st = errno;
  fclose(in);
  gr_free(G);
  errno = st;
  return NULL;
}



/**
   \brief Read an undirected unweighted graph in dimacs format.

   Vertex indices are remapped from [1..n] to [0..n-1].
**/
graph* gr_read_dimacs_uu(char* filename) {

  int i,j,k;

  graph* G = 0;

  FILE* f = fopen(filename,"r");
  if (!f) return 0;

  char* buff = malloc(256*sizeof(char));
  if (!buff) goto ENOMEMH;
  size_t buffs = 256;

  // Gets number of vertices from problem line:
  buff[0] = 0;
  while (buff[0] != 'p' || buff[1] != ' ')
    if (getline(&buff,&buffs,f) == -1) goto EILSEQH;

  k = 2;
  while (buff[k] != ' ')
    k++;

  int st = sscanf(buff+k,"%d %d",&i,&j);
  if (st != 2) goto EILSEQH;

  G = gr_alloc('u',i);
  if (!G) goto ENOMEMH;

  // Gets edges:
  while (getline(&buff,&buffs,f) != -1) {
    if (buff[0] == 'e' && buff[1] == ' ') {
      int st = sscanf(buff+2,"%d %d",&i,&j);
      if (st != 2) goto EILSEQH;
      gr_add_edge(G,i-1,j-1,0);
    }
  }

  fclose(f);
  return G;

 ENOMEMH: k = errno; goto STDH;
 EILSEQH: k = EILSEQ; goto STDH;
 STDH:
  gr_free(G);
  errno = k;
  return 0;
}




/**
  \brief Print a graph to the screen as a list of vertices and edges.
**/
int gr_print(graph *G) {

  printf("Graph type %c\n", G->type);

  printf("Vertices: %d\n", G->n);

  int i;
  for (i=0; i<G->n; i++) {

    printf("vtx %d, flag = %d, N = ",i,G->V[i].flag);

    edge* e = G->V[i].N;
    while (e) {
      printf("%d ",e->term);
      e = e->next;
    }
    printf("\n");
  }

  return 1;
}




/**
  \brief Write a graph to a file in GraphViz dot format.

  Vertices with a non-zero flag field are colored gray.  The others are colored white.

  \param G The graph.
  \param labels An array of l strings that will be used as labels for nodes
   {0,...,l-1} of G.  Nodes {l,l+1,...,n-1} will be labeled l,l+1,...,n-1.  It
   may be NULL.
  \param l The number of labels.
  \param print_edge_weights If non-zero then edge weights are printed.
  \param title A title for the diagram.  May be NULL.
  \param file The output file name.

  \return On success this function returns 1.  On failure it returns 0 and errno
  remains as set by fopen().
**/
int gr_write_dot(graph* G,
                 char** labels, int l,
                 int print_edge_weights, char* title, char* file) {

  char** L = malloc(G->n*sizeof(char*));
  if (!L) { return 0; }

  if (!labels)
    l = 0;

  int i,j;

  for (i=0; i<l; i++)
    L[i] = labels[i];

  for (i=l; i<G->n; i++) {
    if (asprintf(&(L[i]), "%d", i) == -1) {
      for (j=0; j<i; j++)
        free(L[j]);
      free(L);
      return 0;
    }
  }

  int maxlen = 0;
  for (i=0; i<G->n; i++) {
    int len = strlen(L[i]);
    if (len > maxlen)
      maxlen = len;
  }


  FILE *f = fopen(file,"wt");
  if (!f) return 0;


  char conn[3];

  if (G->type == 'd') {
    strcpy(conn,"->");
    fprintf(f,"graph {\n");
  }
  else {
    strcpy(conn,"--");
    fprintf(f,"graph {\n");
  }

  fprintf(f,"node [ fontsize=11 ]\n");
  fprintf(f,"edge [ fontsize=10 ]\n");

  //fprintf(f,"layout=fdp\nstart=1\nk=1.0\nmaxiter=5000\n");
  fprintf(f,"layout=neato\n");

  // From dot manual:
  // dot − filter for drawing directed graphs
  // neato − filter for drawing undirected graphs
  // twopi − filter for radial layouts of graphs
  // circo − filter for circular layout of graphs
  // fdp − filter for drawing undirected graphs
  // sfdp − filter for drawing large undirected graphs
  // patchwork − filter for squarified tree maps
  // osage − filter for array-based layouts

  char shape[10];
  if (maxlen <= 3)
    strcpy(shape,"circle");
  else
    strcpy(shape,"box");

  for (i=0; i<G->n; i++)
    if (G->V[i].flag)
      fprintf(f,"%d [label=\"%s\",style=filled,shape=%s,width=0.1,height=0.1,margin=0.02,fillcolor=gray];\n",i,L[i],shape);
    else
      fprintf(f,"%d [label=\"%s\",shape=%s,width=0.1,height=0.1,margin=0.02];\n",i,L[i],shape);

  for (i=0; i<G->n; i++) {
    edge* e = G->V[i].N;

    if (print_edge_weights) {
      while (e) {
        if (i < e->term)
          fprintf(f,"%d %s %d [label=\"%.3lf\"];\n",i,conn,e->term,e->w);
        e = e->next;
      }
    }
    else {
      while (e) {
        if (i < e->term)
          fprintf(f,"%d %s %d;\n",i,conn,e->term);
        e = e->next;
      }
    }
  }

  if (title) {
    fprintf(f,"labelloc=\"b\"\n");
    fprintf(f,"label=\"%s\"\n",title);
  }

  fprintf(f,"}\n");

  for (i=l; i<G->n; i++)
    free(L[i]);
  free(L);

  fclose(f);

  return 1;
}



/**
  \brief Write a graph to a file in json format.

  The vertices whose flag field is non-zero are colored gray, otherwise white.

  \param G The graph.
  \param print_edge_weights If non-zero then edge weights are printed.
  \param file The output file name.

  \return On success this function returns 1.  On failure it returns 0 and errno
  remains as set by fopen().
**/
int gr_write_json(graph* G, int print_edge_weights, char* file) {

  FILE *f = fopen(file,"wt");
  if (!f) return 0;

  fprintf(f,"var data = {\n");
  fprintf(f,"\"nodes\": [\n");

  int i;
  for (i=0; i<G->n; i++)
    if (G->V[i].flag)
      fprintf(f,"{\"node\": \"L\",\"size\": 2},\n");
    else
      fprintf(f,"{\"node\": \"I\",\"size\": 2},\n");

  fprintf(f,"],\n");
  fprintf(f,"\"links\":[\n");

  for (i=0; i<G->n; i++) {
    edge* e = G->V[i].N;
    while (e) {
      if (i < e->term)
        fprintf(f,"{\"source\": %d, \"target\": %d, \"bond\": 1},\n",i,e->term);
      e = e->next;
    }
  }

  fprintf(f,"]\n");
  fprintf(f,"}\n");

  fclose(f);

  return 1;
}



/**
  \brief Write a graph to a file in GraphViz dot format having image nodes.

  \param G The graph.
  \param image_files Image files for nodes [0,l-1].
  \param l The number of nodes that will be displayed as images.
  \param print_edge_weights If non-zero then edge weights are printed.
  \param file The output file name.

  \return On success this function returns 1.  On failure it returns 0 and errno
  remains as set by fopen().
**/
int gr_write_dot_images(graph* G, char** image_files, int l, int print_edge_weights, char* file) {

  FILE *f = fopen(file,"wt");
  if (!f) return 0;

  fprintf(f,"graph {\n");
  fprintf(f,"node [ shape=circle width=0.2 fillcolor=gray style=filled label=\"\" ]\n");
  fprintf(f,"edge [ label=\"\" ]\n");

  fprintf(f,"layout=fdp\nstart=51\nk=0.75\nmaxiter=5000\n");
  //fprintf(f,"layout=twopi\n");

  int i;
  for (i=0; i<G->n; i++)
    if (i < l)
        fprintf(f,"%d [shape=none,image=\"%s\"];\n",i,image_files[i]);
    else
      fprintf(f,"%d [fixedsize=true];\n",i);


  for (i=0; i<G->n; i++) {
    edge* e = G->V[i].N;

    if (print_edge_weights) {
      while (e) {
        if (i < e->term)
          fprintf(f,"%d -- %d ;\n",i,e->term);
          //fprintf(f,"%d -- %d [label=\"%.4lf\"];\n",i,e->term,e->w);
        e = e->next;
      }
    }
    else {
      while (e) {
        if (i < e->term)
          fprintf(f,"%d -- %d;\n",i,e->term);
        e = e->next;
      }
    }
  }
  fprintf(f,"}\n");

  fclose(f);

  return 1;
}




/**
   \brief Create a Newick representation of a tree.

   \param T The tree.

   \return A string with the Newick representation of T or NULL if memory
   allocation fails somewhere. errno will remain as set by the failing function.
**/
char* tree_newick(graph *T) {

  int i;
  char* r;

  if (T->n == 0)
    return asprintf(&r,";") < 0 ? 0 : r;

  if (T->n == 1)
    return asprintf(&r,"(0);") < 0 ? 0 : r;

  if (T->n == 2)
    return asprintf(&r,"(0,1):%.16e;",gr_edge_weight(T,0,1)) < 0 ? 0 : r;

  int marked[T->n];
  for (i=0; i<T->n; i++)
    marked[i] = 0;


  /* The recusive construction function: */
  char* newick(int k) {

    marked[k] = 1;

    char* r = 0;
    char* s = 0;
    char* t = 0;

    // Count the unmarked neighbors of k:
    int n = 0;
    edge* e = T->V[k].N;
    while (e) {
      if (!marked[e->term])
        n++;
      e = e->next;
    }

    if (n == 0) {
      // All neighbors are marked, k is a trivial tree:
      if (asprintf(&r,"%d",k) < 0) return 0;
    }
    else {
      // Join neigboring trees by k into string r. r starts with a dummy char in it:
      r = calloc(1,sizeof(char));
      if (!r) return 0;

      e = T->V[k].N;
      while (e) {
        if (!marked[e->term]) {
          t = r;
          s = newick(e->term);
          if (!s) {
            free(t);
            return 0;
          }

          r = 0;
          //if (asprintf(&r,"%s,%s:%.16e",t,s,e->w) < 0) {
          if (asprintf(&r,"%s,%s:%.1lf",t,s,e->w) < 0) {
            free(s);
            free(t);
            return 0;
          }
          free(s);
          free(t);
        }
        e = e->next;

        printf("%s\n",r);
      }

      t = r;
      if (asprintf(&r,"(%s)%d",t+1,k) < 0) {
        free(t);
        return 0;
      }
      free(t);
    }

    return r;
  }

  int c1,c2,dia;
  tree_center(T,&c1,&c2,&dia);

  char* t = newick(c1);
  if (!t) return 0;

  if (asprintf(&r,"%s:0.0;",t) < 0) return 0;
  free(t);

  return r;
}




/**
   \brief Create a Newick representation of a tree with given labels.

   \param T The tree.

   \param labels An array of l strings that will be used as labels for nodes
   {0,...,l-1} of T.  Nodes {l,l+1,...,n-1} will be labeled l,l+1,...,n-1.  It may
   be NULL, which implies l=0.

   \param l The number of labels.

   \return A string with the Newick representation of T or NULL if memory
   allocation fails somewhere. errno will remain as set by the failing function.
**/
char* tree_newicks(graph *T, char** labels, int l) {

  int i,j;
  char* r;


  // Define node labels:
  char** L = malloc(T->n*sizeof(char*));
  if (!L) { return 0; }

  if (!labels)
    l = 0;

  // Copy labels for nodes <l:
  for (i=0; i<l; i++) {
    if (asprintf(&(L[i]), "%s", labels[i]) < 0) {
      for (j=0; j<i; j++)
        free(L[j]);
      free(L);
      return 0;
    }
  }

  // Set numeric labels for nodes >=l, <n:
  for (i=l; i<T->n; i++) {
    if (asprintf(&(L[i]), "%d", i) < 0) {
      for (j=0; j<i; j++)
        free(L[j]);
      free(L);
      return 0;
    }
  }


  if (T->n == 0)
    return asprintf(&r,";") < 0 ? 0 : r;

  if (T->n == 1) {
    if (l > 0)
      return asprintf(&r,"(%s)",L[0]) < 0 ? 0 : r;
    else
      return asprintf(&r,"(0)") < 0 ? 0 : r;
  }

  if (T->n == 2) {
    if (l >= 2)
      return asprintf(&r,"(%s,%s):%.16e;",L[0],L[1],gr_edge_weight(T,0,1)) < 0 ? 0 : r;
    else if (l == 1)
      return asprintf(&r,"(%s,1):%.16e;",L[0],gr_edge_weight(T,0,1)) < 0 ? 0 : r;
    else
      return asprintf(&r,"(0,1):%.16e;",gr_edge_weight(T,0,1)) < 0 ? 0 : r;
  }

  int *marked = calloc(T->n,sizeof(int));


  /* The recusive construction function: */
  char* newicks(int k) {

    marked[k] = 1;

    char* r = 0;
    char* s = 0;
    char* t = 0;

    int n = 0;
    edge* e = T->V[k].N;
    while (e) {
      if (!marked[e->term])
        n++;
      e = e->next;
    }

    if (n == 0) {
      // All neighbors are marked, k is a trivial tree:
      if (k < l) {
        if (asprintf(&r,"%s",L[k]) < 0) return 0;
      }
      else {
        if (asprintf(&r,"%d",k) < 0) return 0;
      }
    }
    else {
      // Join neigboring trees by k into string r. r starts with a dummy char in it:
      r = calloc(1,sizeof(char));
      if (!r) return 0;

      e = T->V[k].N;
      while (e) {
        if (!marked[e->term]) {
          t = r;
          s = newicks(e->term);
          if (!s) {
            free(t);
            return 0;
          }

          r = 0;
          if (asprintf(&r,"%s,%s:%.16e",t,s,e->w) < 0) {
            free(s);
            free(t);
            return 0;
          }

          free(s);
          free(t);
        }
        e = e->next;
      }

      t = r;
      if (k < l) {
        if (asprintf(&r,"(%s)%s",t+1,L[k]) < 0) {
          free(t);
          return 0;
        }
      }
      else {
        //if (asprintf(&r,"(%s)%d",t+1,k) < 0) exit(0);
        if (asprintf(&r,"(%s)",t+1) < 0) {
          free(t);
          return 0;
        }
      }

      free(t);
    }

    return r;
  }

  int c1,c2,dia;
  tree_center(T,&c1,&c2,&dia);

  char* t = newicks(c1);
  if (asprintf(&r,"%s:0.00;",t) < 0) return 0;
  free(t);

  free(marked);
  for (i=0; i<T->n; i++)
    free(L[i]);
  free(L);

  return r;
}




/**
   \brief Write the Newick representation of a tree to a file.

   \param T The tree.
   \param labels An array of l strings that will be used as labels for nodes
   {0,...,l-1} of T.  Nodes {l,l+1,...,n-1} will be labeled l,l+1,...,n-1.  It may
   be NULL, which implies l=0.
   \param l The number of labels.
   \param file The output file name.

   \return On success 1 or 0 if a memory allocation or file operation fails.
   errno will remain as set by the failing function.
**/
int tree_write_newicks(graph *T, char** labels, int l, char* file) {

  char* newick = tree_newicks(T,labels,l);
  if (!newick) return 0;

  FILE *f = fopen(file,"wt");
  if (!f) {
    free(newick);
    return 0;
  }
  fprintf(f,"%s\n",newick);
  fclose(f);

  free(newick);

  return 1;
}



/**
   \brief Find the vertex with maximum degree (out-degree) and smallest index
   in G.

   This function does not modify graph flags.

   \param G The graph.
   \param maxdeg The maximum degree (out-degree) in the graph.
   \param v The vertex with maximum degree and smallest index.
**/
void gr_max_degree(graph* G, int* maxdeg, int* v) {

  *maxdeg = 0;

  int i;
  for (i=0; i<G->n; i++) {

    int deg = 0;
    edge* e = G->V[i].N;
    while (e) {
      deg++;
      e = e->next;
    }

    if (deg > *maxdeg) {
      *v = i;
      *maxdeg = deg;
    }
  }
}




/**
   \brief Evaluate the degree (out-degree) of each vertex v.

   Evaluate the degree (non-oriented) or out-degree (oriented) of every vertex v
   and store in v.flag.

   \param G The graph.
   \param maxdeg The maximum degree (out-degree) in the graph.
   \param v The vertex with maximum degree and smallest index.
**/
void gr_eval_degrees(graph* G, int* maxdeg, int* v) {

  *maxdeg = 0;

  int i;
  for (i=0; i<G->n; i++) {

    G->V[i].flag = 0;
    edge* e = G->V[i].N;
    while (e) {
      G->V[i].flag++;
      e = e->next;
    }

    if (G->V[i].flag > *maxdeg) {
      *v = i;
      *maxdeg = G->V[i].flag;
    }
  }
}





//// A datatype to store a reference to edge and its weigh, and a
//// function to sort them.  Needed by MST.
//
//typedef struct redge {
//  int u;    // Edge starting end.
//  edge *e;
//} redge;
//
//
//int compare(const void *e1, const void *e2) {
//
//  redge *e = (redge*) e1;
//  redge *f = (redge*) e2;
//
//  if (e->e->w < f->e->w)
//    return -1;
//  else if (e->e->w > f->e->w)
//    return 1;
//  else
//    if (e->u < f->u)
//      return -1;
//    else if (e->u > f->u)
//      return 1;
//    else
//      if (e->e->v < f->e->v)
//        return -1;
//      else if (e->e->v > f->e->v)
//        return 1;
//      else
//        return 0;
//}
//
//
///*
//Kruskal's MST.
//
//Returns 1 for success, 0 otherwise.  Sets w equal to thetree weight.
//Sets in=1 for tree edges and 0 for non-tree edges.  For connected
//simple graphs sets unique=0/1 for nonunique/unique MST.
//*/
//
//int MST(graph *G, int *w, int *unique) {
//
//  int i,u;
//  edge *a;
//
//  // Builds a list of edges for sorting:
//  redge *V = malloc(G->m*sizeof(redge));
//  if (!V)
//    return 0;
//
//  i=0;
//  for (u=1; u<=G->n; u++) {
//    a = G->V[u].N;
//    while (a) {
//      if (G->type == 0 || a->v > u) {
//        V[i].u = u;
//        V[i++].e = a;
//      }
//      a->in = 0;
//      a = a->prox;
//    }
//  }
//
//  qsort(V, G->m, sizeof(redge), compare);
//
//  // Kruskal with disjoint sets and uniqueness test:
//  disjoint_sets *S = init_ds(G->n);
//  int last;
//
//  // Adds n-2 edges to the tree and records the last:
//  *w = 0;
//  last = 0;
//  i=0;
//  while (!last) {
//    if (find_ds(S,V[i].u) != find_ds(S,V[i].e->v)) {
//      if (S->n == 2)
//        last = i;
//      else {
//        union_ds(S,V[i].u,V[i].e->v);
//        *w += V[i].e->w;
//        V[i].e->in = 1;
//      }
//    }
//    i++;
//  }
//
//  // Tests for valid edges with weight equal to the last:
//  *unique = 1;
//  if (G->m >= G->n) {
//    i = last + 1;
//    while (i<G->m && V[last].e->w == V[i].e->w) {
//      if (find_ds(S,V[i].u) != find_ds(S,V[i].e->v)) {
//        *unique = 0;
//        i = G->m;
//      }
//      i++;
//    }
//  }
//
//  // Adds the last edge:
//  *w += V[last].e->w;
//  V[last].e->in = 1;
//
//  return 0;
//}
//
//











//int* dfs(graph G, int v, int *acyclic) {
//
//  int *p;  // predecessors array.
//  int *color; // array to hold colors. 0 for white, 1 for gray, 2 for black.
//  int *tf, *tl;  // arrays to hold timestamps for vertices.
//
//  stack S;  // An stack for edges. (v,0) is used to indicated that
//          // N(v) is fully visited
//
//  S = init_stack(n);
//
//  push(S,v);
//  push(S,0);
//
//  a = G->V[v].N;
//  while (a) {
//    push(S,v);
//    push(S,a->v);
//    a = a->prox;
//  }
//
//  while (!empty(S)) {
//
//    // Take edge (u,v) from stack:
//    v = push(S);
//    u = push(S);
//
//    if (v == 0) {
//      tl[v] = time++;
//      color[v] = 2;
//    }
//    else {
//      if (color[v] == 0) {
//      p[v] = u;
//      tf[v] = time++;
//      color[v] = 1;
//      }
//      a = G->V[v].N;
//      while (a) {
//      if (color[a->v] == 0) {
//        push(S,v);
//        push(S,a->v);
//        a = a->prox;
//      }
//      }
//    }
//  }
//
//  destroy_stack(S);
//  free(color);
//  free(tf);
//  free(tl);
//}



/**
   \brief A BFS from s that stores the distances in the flag field of each vertex.

   \return On success it returns 1.  Otherwise it returns 0 and errno is left as
   set by malloc() on failure.
**/
int gr_bfs(graph* G, int s) {

  dynarray* Q = da_alloc('i',max(10,G->m/G->n*G->m/G->n),0);
  if (!Q) return 0;

  int i;
  for (i=0; i<G->n; i++)
    G->V[i].flag = INT_MAX;

  da_push(Q,s);
  G->V[s].flag = 0;

  while (da_size(Q)) {
    int u;
    da_shift(Q,&u);

    edge* a = G->V[u].N;
    while (a) {
      if (G->V[a->term].flag == INT_MAX) {
        da_push(Q,a->term);
        G->V[a->term].flag = G->V[u].flag + 1;
      }
      a = a->next;
    }
  }

  da_free(Q);
  return 1;
}





/**
   \brief Evaluate the diameter of a graph.

   \internal Implemented as multiple BFSs.
**/
int gr_diameter(graph* G) {

  dynarray* Q = da_alloc('i',max(10,G->m/G->n*G->m/G->n),0);
  if (!Q) return 0;

  int dia = 0, u, v;

  for (v=0; v<G->n; v++) {

    for (u=0; u<G->n; u++)
      G->V[u].flag = INT_MAX;

    da_push(Q,v);
    G->V[v].flag = 0;

    while (da_size(Q)) {
      int u;
      da_shift(Q,&u);

      edge* a = G->V[u].N;
      while (a) {
        if (G->V[a->term].flag == INT_MAX) {
          da_push(Q,a->term);
          G->V[a->term].flag = G->V[u].flag + 1;
          if (G->V[a->term].flag > dia)
            dia = G->V[a->term].flag;
        }
        a = a->next;
      }
    }
  }

  da_free(Q);
  return dia;
}




/**
   \brief All-pairs shortest path weigths.

   This is an inplementation of the Floyd-Warshal algorithm.

   \return A new n x n array of pairwise shortest path weights, DBL_MAX meaning unreachable.
**/
double** gr_apspw(graph* G) {

  int m,i,j;
  double** M = (double**) m_alloc('d',G->n,G->n);
  if (!M) return 0;

  for (i=0; i<G->n; i++) {
    for (j=0; j<G->n; j++)
      M[i][j] = INFINITY;
    M[i][i] = 0;
  }

  for (i=0; i<G->n; i++) {
    edge* e = G->V[i].N;
    while (e) {
      M[i][e->term] = e->w;
      e = e->next;
    }
  }

  for (m=0; m<G->n; m++)
    for (i=0; i<G->n; i++)
      for (j=0; j<G->n; j++)
        if (isfinite(M[i][m]) && isfinite(M[m][j]) && M[i][m] + M[m][j] < M[i][j])
          M[i][j] = M[i][m] + M[m][j];
          // P[i][j] = P[m][j];

  return M;
}




/**
   \brief Evaluate all-pairs distances among vertices {0,...,k-1}.

   Implemented as multiple BFSs.

   \return A new k x k array with pairwise distances, INT_MAX meaning unreachable.
**/
int** gr_apd(graph* G, int k) {

  dynarray* Q;
  int** M = 0;
  int* mark = 0;

  Q = da_alloc('i',max(10,G->m/G->n*G->m/G->n),0);
  if (!Q) return 0;

  M = (int**) m_alloc('i',k,k);
  if (!M) goto nomem;

  mark = malloc(G->n*sizeof(int));
  if (!mark) goto nomem;

  int i,j,s;

  for (i=0; i<k; i++) {
    for (j=0; j<k; j++)
      M[i][j] = INT_MAX;
    M[i][i] = 0;
  }

  for (s=0; s<k; s++) {

    for (i=0; i<G->n; i++)
      mark[i] = INT_MAX;

    da_push(Q,s);
    mark[s] = 0;
    M[s][s] = 0;

    while (da_size(Q)) {
      int u;
      da_shift(Q,&u);

      edge* a = G->V[u].N;
      while (a) {
        if (mark[a->term] == INT_MAX) {
          da_push(Q,a->term);
          mark[a->term] = mark[u] + 1;
          if (a->term < k)
            M[s][a->term] = mark[a->term];
        }
        a = a->next;
      }
    }
  }

  da_free(Q);
  free(mark);
  return M;

 nomem:
  da_free(Q);
  m_free((void**)M,k);
  free(mark);

  return 0;
}




/**
   \brief Single-source paths and path weights in a tree.

   This function will change the contents of flag in each vertex.

   \param G An acyclic graph.
   \param s The source.
   \param pi The predecessor of each vertex in the path to s.  The predecessor
   of s will be -1. Unreachable vertices' predecessors will be -1.
   \param weights The path weights. Unreachable vertices' weight will be DBL_MAX.

   \return On success it returns 1.  Otherwise it returns 0 and errno is left as
   set by malloc() on failure.
**/
int tree_sspw(graph *T, int s, int* pi, double* weights) {

  dynarray* Q = da_alloc('i',max(10,T->m/T->n*T->m/T->n),0);
  if (!Q) return 0;

  int i;
  for (i=0; i<T->n; i++) {
    T->V[i].flag = 0;
    pi[i] = -1;
    weights[i] = DBL_MAX;
  }

  da_push(Q,s);
  T->V[s].flag = 1;
  weights[s] = 0;
  pi[s] = -1;

  while (da_size(Q)) {
    int u;
    da_shift(Q,&u);

    edge* a = T->V[u].N;
    while (a) {
      if (T->V[a->term].flag == 0) {
        da_push(Q,a->term);
        weights[a->term] = weights[u]+a->w;
        pi[a->term] = u;
        T->V[a->term].flag = 1;
      }
      a = a->next;
    }
  }

  da_free(Q);
  return 1;
}






/**
   \brief Find the maximal cliques in an undirected graph.

   Uses Bron-Kerbosch maximal cliques enumeration (algorithm 1, bactraking).

   \return A list of pointers to cliques, each clique an array of int whose first
   position holds its length.
**/
list* gr_max_cliques1(graph* G) {

  int i;
  list* cliques = list_alloc('v');

  int v, maxoutdeg;
  gr_eval_degrees(G,&maxoutdeg,&v);

  // A stack with its first free position in stk[0]:
  int* stk = malloc((maxoutdeg+1)*sizeof(int));
  if (!stk) die("gr_max_cliques1 stk");
  stk[0] = 1;


  /**
     Bron-Kerbosch maximal cliques enumeration (algorithm 1, bactraking).

     The compsub set is stored in the stack stk.  Not and candidates sets are
     implemented in the array NC.  Position 0 of NC is not used (it means empty),
     index ne points to the last element in not and index ce points to the last
     element in candidates.

     Adds cliques to the list cliques.

     C. Bron and J. Kerbosch. Algorithm 457: Finding all cliques of an
     undirected graph Comm. of the ACM, 16, 575-577,1973.
  **/
  void BK1(int* NC, int ne, int ce) {

    int c,i,j;

    for (c=ne+1; c<=ce; c++) {
      // Selects a candidate: just get c.

      // Adds the candidate to the growing clique:
      stk[stk[0]++] = NC[c];

      // Creates new not and candidates sets:
      int* newNC = malloc((G->V[NC[c]].flag+1)*sizeof(int));
      if (!newNC) die("BK2 newNC");

      // N' = N \cap N[u]:
      int newne = 0;
      for (j=1; j<=ne; j++)
        if (gr_edge(G,NC[j],NC[c]))
          newNC[++newne] = NC[j];

      // C' = C \cap N[u]:
      int newce = newne;
      for (j=ne+2; j<=ce; j++)
        if (gr_edge(G,NC[j],NC[c]))
          newNC[++newce] = NC[j];

      if (newce == 0) { // A clique was found.
        int* C = malloc((stk[0])*sizeof(int));
        if (!C) die("BK2 C");

        C[0] = stk[0]-1;
        for (i=1; i<stk[0]; i++)
          C[i] = stk[i];
        list_push(cliques,1,(void*)C);
      }
      else if (newne < newce) // Any candidate is left.
        BK1(newNC,newne,newce);

      // Removes the candidate from the growing clique:
      stk[0]--;

      // Adds to not:
      ne++;

      free(newNC);
    }
  }

  int* NC = malloc((G->n+1)*sizeof(int));
  if (!NC) die("gr_max_cliques1 NC");

  NC[0] = 0;
  for (i=0; i<G->n; i++)
    NC[i+1] = i;

  BK1(NC,0,G->n);

  free(NC);
  free(stk);

  return cliques;
}




/**
   \brief Find the maximal cliques in an undirected graph.

   Uses Bron-Kerbosch maximal cliques enumeration (algorithm 2, bound with
   disconnections counting).

   \return A list of pointers to cliques, each clique an array of int whose first
   position holds its length.
**/
list* gr_max_cliques2(graph* G) {

  int i;
  list* cliques = list_alloc('v');

  int v,maxoutdeg;
  gr_eval_degrees(G,&maxoutdeg,&v);

  // A stack with its first free position in stk[0]:
  int* stk = malloc((maxoutdeg+1)*sizeof(int));
  if (!stk) die("gr_max_cliques2 stk");
  stk[0] = 1;


  /**
     Bron-Kerbosch maximal cliques enumeration (algorithm 2, bound with
     disconnections counting).

     The compsub set is stored in the stack stk. Not and candidates sets are
     implemented in the array NC. Position 0 of NC is not used, index ne points
     to the last element in not and index ce points to the last element in
     candidates.

     Adds cliques to the list cliques.

     C. Bron and J. Kerbosch. Algorithm 457: finding all cliques of an
     undirected graph Comm. of the ACM, 16, 575-577,1973.
  **/
  void BK2(int* NC, int ne, int ce) {

    int i, j;

    // Finds a candidate with the minimum number of disconections and the
    // fixed-point:
    int fixp = 0;
    int s = 0;
    int nod = 0;
    int minnod = ce;
    i = 1;

    while (i <= ce && minnod != 0) {
      int count = 0, pos = 0;
      j = ne + 1;
      while (j <= ce && count < minnod) {
        if (i == j)
          ;
        else {
          if (!gr_edge(G,NC[i],NC[j])) {
            count++;
            pos = j;
          }
        }
        j++;
      }

      if (count < minnod) {
        fixp = NC[i];
        minnod = count;
        // If fixed point is initially selected from C, it will be the first to
        // be added to not and then nod will be increased by 1. Otherwise if the
        // fixed point is in N, then the candidate is the last disconnection of
        // the fixed point in C.
        if (i <= ne)
          s = pos;
        else {
          s = i;
          nod = 1;
        }
      }

      i++;
    }

    nod += minnod;
    while (nod >= 1) {
      // Interchange c with ne+1:
      int aux = NC[s];
      NC[s] = NC[ne + 1];
      NC[ne + 1] = aux;

      s = ne + 1;

      // Creates new set not:
      int* newNC = malloc((G->V[NC[s]].flag+1)*sizeof(int));
      if (!newNC) die("BK2 newNC");

      // N' = N \cap N[u]:
      int newne = 0;
      for (i = 1; i <= ne; i++)
        if (gr_edge(G,NC[i],NC[s]))
          newNC[++newne] = NC[i];

      // C' = C \cap N[u]:
      int newce = newne;
      for (i = ne + 2; i <= ce; i++)
        if (gr_edge(G,NC[i],NC[s]))
          newNC[++newce] = NC[i];

      // Adds the candidate to the growing clique:
      stk[stk[0]++] = NC[s];

      // A maximal clique was found:
      if (newce == 0) {
        int* C = malloc((stk[0])*sizeof(int));
        if (!C) die("BK2 C");
        C[0] = stk[0]-1;
        for (i=1; i<stk[0]; i++)
          C[i] = stk[i];
        list_push(cliques,1,(void*)C);
      }
      else if (newne < newce)
        BK2(newNC,newne,newce);

      // Removes the candidate from the growing clique:
      stk[0]--;

      // Add to not:
      ne++;

      free(newNC);

      // Selects a candidate disconnected to the fixed point:
      if (nod > 1) {
        s = ne + 1;
        while (gr_edge(G,fixp,NC[s]))
          s++;
      }
      nod--;
    }
  }

  int* NC = malloc((G->n+1)*sizeof(int));
  if (!NC) die("gr_max_cliques2 NC");

  NC[0] = 0;
  for (i=0; i<G->n; i++)
    NC[i+1] = i;

  BK2(NC,0,G->n);

  free(NC);
  free(stk);

  return cliques;
}



/**
   \brief Find centers and the diameter of a tree.

   This function will change the contents of flag in each vertex.

   \param G An acyclic undirected graph.
   \param c1 The index of the center.
   \param c2 The index of the other center.  If the center is unique, c2=c1.
   \param The diameter.

   \return On success it returns 1.  Otherwise it returns 0 and errno is left as
   set by malloc() on failure.
**/
int tree_center(graph *T, int* c1, int* c2, int *dia) {

  int i;
  int s=0, t=0;

  // BFS from 0:
  dynarray* Q = da_alloc('i',max(10,T->m/T->n*T->m/T->n),0);
  if (!Q) return 0;

  for (i=0; i<T->n; i++)
    T->V[i].flag = INT_MAX;

  da_push(Q,0);
  T->V[0].flag = 0;

  while (da_size(Q)) {
    int u;
    da_shift(Q,&u);
    s = u;

    edge* a = T->V[u].N;
    while (a) {
      if (T->V[a->term].flag == INT_MAX) {
        da_push(Q,a->term);
        T->V[a->term].flag = T->V[u].flag + 1;
      }
      a = a->next;
    }
  }

  // BFS from s and record paths:
  int* pi = malloc(T->n*sizeof(int));
  if (!pi) {
    da_free(Q);
    return 0;
  }

  for (i=0; i<T->n; i++)
    T->V[i].flag = INT_MAX;

  da_push(Q,s);
  T->V[s].flag = 0;
  pi[s] = -1;

  while (da_size(Q)) {
    int u;
    da_shift(Q,&u);
    t = u;

    edge* a = T->V[u].N;
    while (a) {
      if (T->V[a->term].flag == INT_MAX) {
        da_push(Q,a->term);
        T->V[a->term].flag = T->V[u].flag + 1;
        pi[a->term] = u;
      }
      a = a->next;
    }
  }

  *dia = T->V[t].flag;

  i = t;
  while (T->V[i].flag > *dia/2+1)
    i = pi[i];

  if (*dia % 2 == 1) {
    *c1 = i;
    *c2 = pi[i];
  }
  else
    *c1 = *c2 = pi[i];

  free(pi);
  da_free(Q);

  return 1;
}



/**
   \brief Find a leaf that is nearest to a given node of a tree.

   This function will not change the value of any flag fields.

   \param T An acyclic graph.
   \param s The source.
   \param l A nearest leaf.
   \param d The distance to l.

   \return On success it returns 1.  Otherwise it returns 0 and errno is left as
   set by malloc() on failure.
**/
int tree_nearest_leaf(graph *G, int s, int* l, int* d) {

  dynarray* Q = da_alloc('i',max(10,G->m/G->n*G->m/G->n),0);
  if (!Q) return 0;

  int* mark = malloc(G->n*sizeof(int));
  if (!mark) {
    da_free(Q);
    return 0;
  }

  int i;
  for (i=0; i<G->n; i++)
    mark[i] = -1;

  da_push(Q,s);
  mark[s] = 0;

  while (da_size(Q)) {
    da_shift(Q,&i);

    int n = 0;
    edge* a = G->V[i].N;
    while (a) {
      if (mark[a->term] == -1) {
        da_push(Q,a->term);
        mark[a->term] = mark[i]+1;
      }
      n++;
      a = a->next;
    }

    if (n == 1) {
      *l = i;
      *d = mark[i];
      da_free(Q);
      free(mark);
      return 1;
    }
  }

  da_free(Q);
  free(mark);

  return 1;
}




/**
   \brief All-pairs path weights among vertices [0,k-1].

   Implemented as multiple BFSs.  It doesn't modify the flag fields in vertices.

   \return A new k x k ragged array with pairwise distances, DBL_MAX meaning
   unreachable.  On error it returns NULL and errno remains as set by malloc.
**/
double** tree_apw(graph* G, int k) {

  dynarray* Q;
  double** M = 0;
  double* W = 0;

  Q = da_alloc('i',G->n,0);
  if (!Q) return 0;

  M = (double**) m_alloc('d',k,k);
  if (!M) goto nomem;

  W = malloc(G->n*sizeof(double));
  if (!W) goto nomem;

  int i,j,s;

  for (i=0; i<k; i++)
    for (j=0; j<k; j++)
      M[i][j] = DBL_MAX;

  for (s=0; s<k; s++) {
    for (i=0; i<G->n; i++)
      W[i] = DBL_MAX;

    da_push(Q,s);
    W[s] = 0.0;
    M[s][s] = 0.0;

    while (da_size(Q)) {
      int u;
      da_shift(Q,&u);

      edge* a = G->V[u].N;
      while (a) {
        if (W[a->term] == DBL_MAX) {
          da_push(Q,a->term);
          W[a->term] = W[u] + a->w;
          if (a->term < k)
            M[s][a->term] = W[a->term];
        }
        a = a->next;
      }
    }
  }

  da_free(Q);
  free(W);
  return M;

 nomem:
  i = errno;
  da_free(Q);
  free(W);
  m_free((void**)M,k);
  errno = i;
  return 0;
}


/**
   \brief Get edges of an undirected graph whose both ends have degree larger
   that 1.

   \param G The graph.

   \param both_directions If both_directions is 1 then for each edge (u,v) both
   u,v and v,u are returned, and if both_directions is 0 then for each edge
   (u,v) only u,v is returned.

   \param edges A new array with 2*m vertex indices, such that
   (edges[0], edges[1]), ..., (edges[2m-2], edges[2m-1])
   are the internal edges.

   \param m The number of internal edges.

   \return On success it returns 1 and edges points to a new array with 2*m
   elements.  On failure it returns 0, m gets 0, and errno remains as set by
   malloc() on failure.
**/
int gr_internal_edges(graph *G, int both_directions, int** edges, int* m) {

  int i, n;

  gr_eval_degrees(G, &i, &n);

  n = 0;
  for (i=0; i<G->n; i++) {
    if (G->V[i].flag > 1) {

      edge* e = G->V[i].N;
      while (e) {
        if ((both_directions || i < e->term) && G->V[e->term].flag > 1)
          n++;
        e = e->next;
      }
    }
  }

  int* E = malloc(2*n*sizeof(int));
  if (!edges) {
    *m = 0;
    return 0;
  }

  *m = n;

  n = 0;
  for (i=0; i<G->n; i++) {

    if (G->V[i].flag > 1) {

      edge* e = G->V[i].N;
      while (e) {
        if ((both_directions || i < e->term) && G->V[e->term].flag > 1) {
          E[n++] = i;
          E[n++] = e->term;
        }
        e = e->next;
      }
    }
  }

  *edges = E;
  return 1;
}
