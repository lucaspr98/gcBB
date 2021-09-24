
void bwsd(char* file1, char* file2, size_t n, int k, double *expectation, double *entropy, int mem, int printBoss, char coverage_type, size_t total_coverage, short complement);

size_t apply_coverage_merge(int primaryCoverage, int secondaryCoverage, size_t *rl_freq, size_t pos);

double bwsd_expectation(size_t *t, size_t s, size_t n);

double bwsd_shannon_entropy(size_t *t, size_t s, size_t n);

void print_bwsd_matrixes(double **Dm, double **De, char **files, int files_n, char *path, int k, char coverage_type);

double log2(double i);