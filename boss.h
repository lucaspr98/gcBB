
void Wi_sort(char *Wi, int *Wm, int *colors, int *coverage, int start, int end);

int boss_construction(int *LCP, int *DA, char *BWT, int *C, int *last, char *W, int *Wm, int *colors, int n, int k, int samples, int *reduced_LCP, int *coverage, int *total_coverage);

void print_boss_result(int boss_len, int id1, int id2, char *file1, char *file2, int *C, int *last, char *W, int *Wm, int *colors, int *reduced_LCP, int *coverage, int total_coverage);