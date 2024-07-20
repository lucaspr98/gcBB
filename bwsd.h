
void bwsd(char* file1, char* file2, int k, double *expectation, double *entropy, int mem, int printBoss, int consider1, int consider2);

void bwsdAll(char* path, int samples, int k, int mem, double** Dm, double** De);

void applyCoverageMerge(int zeroCoverage, int oneCoverage, size_t *rlFreq, size_t *pos);

double bwsdExpectation(size_t *t, size_t s, size_t n);

double bwsdShannonEntropy(size_t *t, size_t s, size_t n);

double log2(double i);