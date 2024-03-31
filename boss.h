void WiSort(char *Wi, short *Wm, short *colors, int *coverage, short *summarizedSL, int start, int end);

void fixWiLCP(char *W, short *summarizedLCP, int k, int WiSize);

void bossConstruction(FILE *mergeLCP, FILE *mergeDA, FILE *mergeBWT, FILE *mergeSL, size_t n, int k, int samples, int mem, char* file1, char* file2, int printBoss);

/* edgeStatus:
   0: any outgoing edge besides the last one
   1: just one outgoing edge
   2: last outgoing edge from a set
 */
void addEdge(char *W, short **last, short *colors, short *summarizedLCP, short *summarizedSL, int freq, short *Wm, char bwt, char da, short lcp, short sl, int WiSize, int edgeStatus);

