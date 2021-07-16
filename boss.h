
void Wi_sort(char *Wi, int *Wm, int *colors, int *coverage, int start, int end);

int boss_construction(short *LCP, int *DA, char *BWT, int *C, int *last, char *W, int *Wm, int *colors, int n, int k, int samples, short *reduced_LCP, int *coverage, int *total_coverage);

void print_boss_result(int boss_len, char *file1, char *file2, int *C, int *last, char *W, int *Wm, int *colors, short *reduced_LCP, int *coverage, int total_coverage);

/* edge_status:
   0: any outgoing edge besides the last one
   1: just one outgoing edge
   2: last outgoing edge from a set
 */
void add_edge(int i, char *W, int **last, int *colors, short *reduced_LCP, int freq, int *Wm, char bwt, int da, short lcp, int Wi_size, int edge_status);