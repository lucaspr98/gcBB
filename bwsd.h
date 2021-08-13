
void bwsd(char* file1, char* file2, size_t n, int k, double *expectation, double *entropy, int mem, int printBoss);

size_t apply_coverage(short primaryColor, short secondaryColor, int primaryCoverage, int secondaryCoverage, short *rl_color, int *rl_freq, size_t pos);

double bwsd_expectation(int *t, int s, int n);

double bwsd_shannon_entropy(int *t, int s, int n);

void print_bwsd_matrixes(double **Dm, double **De, char **files, int files_n, char *path, int k);