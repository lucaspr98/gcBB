
void Wi_sort(char *Wi, short *Wm, short *colors, int *coverage, int start, int end);

size_t boss_construction(FILE *mergeLCP, FILE *mergeDA, FILE *mergeBWT, FILE *mergeSL, size_t n, int k, int samples, int mem, char* file1, char* file2, int printBoss, size_t *total_coverage);

/* edge_status:
   0: any outgoing edge besides the last one
   1: just one outgoing edge
   2: last outgoing edge from a set
 */
void add_edge(char *W, short **last, short *colors, short *summarized_LCP, short *summarized_SL, int freq, short *Wm, char bwt, int da, short lcp, short sl, int Wi_size, int edge_status);

