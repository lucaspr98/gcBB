void Wi_sort(char *Wi, short *Wm, short *colors, int *coverage, short *summarized_SL, int start, int end);

void fix_Wi_LCP(char *W, short *summarized_LCP, int k, int Wi_size);

size_t boss_construction(FILE *mergeLCP, FILE *mergeDA, FILE *mergeBWT, FILE *mergeSL, size_t n, int k, int samples, int mem, char* path, int printBoss, size_t *totalSampleCoverageInBoss, size_t *totalSampleColorsInBoss);

/* edge_status:
   0: any outgoing edge besides the last one
   1: just one outgoing edge
   2: last outgoing edge from a set
 */
void add_edge(char *W, short **last, short *colors, short *summarized_LCP, short *summarized_SL, int freq, short *Wm, char bwt, int da, short lcp, short sl, int Wi_size, int edge_status);

