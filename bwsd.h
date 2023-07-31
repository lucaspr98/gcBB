
void bwsd(char* path, size_t n, int k, double *expectation, double *entropy, int mem, size_t total_coverage, int consider1, int consider2);

void bwsd_all(char* path, int samples, size_t n, size_t *sample_size, int k, int mem, double** Dm, double** De);

void apply_coverage_merge(int zeroCoverage, int oneCoverage, size_t *rl_freq, size_t *pos);

double bwsd_expectation(size_t *t, size_t s, size_t n);

double bwsd_shannon_entropy(size_t *t, size_t s, size_t n);

double log2(double i);