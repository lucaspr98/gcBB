
void bwsd(int *colors, short *reduced_LCP, int *coverage, size_t n, int k, double *expectation, double *entropy, int mem);

double bwsd_expectation(int *t, int s, int n);

double bwsd_shannon_entropy(int *t, int s, int n);

void print_bwsd_matrixes(double **Dm, double **De, char **files, int files_n, char *path);