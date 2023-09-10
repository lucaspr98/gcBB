void computeFile(char *path, char *file, int memory);

void computeMergeFileAll(char *path, char **files, int numberOfFiles, int memory);

void computeMergeFiles(char *path, char *file1, char *file2, int memory);

void printDistanceMatrixes(double **Dm, double **De, char **files, int files_n, char *path, int k);

// If ALL_VS_ALL, pass path as file1 and NULL as file2
/* update
    0: creates file "w";
    1: append to file "+1";
*/
FILE* getInfoFile(char* file1, char* file2, int k, int update);