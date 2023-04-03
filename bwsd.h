
void bwsd(char* file1, char* file2, size_t n, int k, double *expectation, double *entropy, int mem, int printBoss, size_t total_coverage);

void apply_coverage_merge(int zeroCoverage, int oneCoverage, size_t *rl_freq, size_t *pos, size_t *max_freq);

double bwsd_expectation(size_t *t, size_t s, size_t n);

double bwsd_shannon_entropy(size_t *t, size_t s, size_t n);

double log2(double i);