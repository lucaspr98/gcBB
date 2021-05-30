
void bwsd(int *DA, int n, double *expectation, double *entropy);

double bwsd_expectation(int *t, int s, int n);

double bwsd_shannon_entropy(int *t, int s, int n);

void print_bwsd_matrixes(double **Dm, double **De, char **files, int files_n);