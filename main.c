#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <dirent.h>
#include <unistd.h>
#include <time.h>
#include <libgen.h>

#include "bwsd.h"
#include "boss.h"
#include "external.h"
#include "lib/rankbv.h"

#define FILE_PATH 1024

#ifndef COVERAGE
	#define COVERAGE 0
#endif

#ifndef ALL_VS_ALL
	#define ALL_VS_ALL 0
#endif

#ifndef DEBUG
	#define DEBUG 0
#endif

int compareFiles(const void *element1, const void *element2) {
    const char **file1 = (const char **)element1;
    const char **file2 = (const char **)element2;
    return strcmp(*file1, *file2);
}

char* getPathDirName(char *path, int len){
    if(path[len-1] == '/'){
        path[len-1] = '\0';
        len--;
    }
    int lastSlashPos = 0;
    int i;
    for(i = len-1; i >= 0; i--){
        if(path[i] == '/'){
            lastSlashPos = i;
            break;
        }
    }
    char* dir = calloc(len+1, sizeof(char));
    int j = 0;
    lastSlashPos = path[lastSlashPos] == '/' ? lastSlashPos + 1 : lastSlashPos;
    for(i = lastSlashPos; i <= len; i++){
        dir[j] = path[i];
        j++;
    }
    free(path);
    return dir;
}

int main(int argc, char *argv[]){
    int i, j;
    char **files = (char**)calloc(512, sizeof(char*));
    int k = 32;
    int numberOfFiles = 0;
    char *path;
    int pathLen;
    int opt;
    int memory = 2048;
    int printBoss = 0;

    /******** Check arguments ********/
    int validOpts = 0;
    while ((opt = getopt (argc, argv, "pk:m:")) != -1){
        switch (opt){
            case 'p':
                validOpts+=1;
                printBoss = 1;
                break;
            case 'k':
                validOpts += 2;
                k = atoi(optarg);
                break;
            case 'm':
                validOpts += 2;
                memory = atoi(optarg);
                break;
            case '?':
                if(opt == 'k')
                    fprintf (stderr, "Option -%c requires a integer value.\n", opt);
                else if(opt == 'm')
                    fprintf (stderr, "Option -%c requires a integer value.\n", opt);
                else if (isprint (opt))
                    fprintf (stderr, "Unknown option `-%c'.\n", opt);
                else
                    fprintf (stderr, "Unknown option character `\\x%x'.\n", opt);
                return 1;
            default:
                abort ();
        }
    }

    if(argc-validOpts == 2){
        DIR *folder;
        struct dirent *entry;
        int len;
        pathLen = strlen(argv[argc-1]);
        path = (char*)malloc((pathLen+1)*sizeof(char));
        strncpy(path, argv[argc-1], pathLen+1);

        folder = opendir(argv[argc-1]);
        if(folder == NULL){
            fprintf (stderr, "Unable to read directory %s\n", argv[argc-1]);
            exit(-1);
        }
        while((entry=readdir(folder)) != NULL){
            char *isFastq = strstr(entry->d_name, ".fastq");
            char *isFasta = strstr(entry->d_name, ".fasta");

            if((isFastq && strlen(isFastq) == 6) || (isFasta && strlen(isFasta) == 6)){
                len = strlen(entry->d_name)+1;
                files[numberOfFiles] = (char*)malloc((pathLen+len+2)*sizeof(char));

                strncpy(files[numberOfFiles], entry->d_name, pathLen+len+2);
                numberOfFiles++;
            }
        }

        closedir(folder);
    } else if(argc-validOpts == 4){
        int fileLen;
        pathLen = strlen(argv[argc-3]);
        path = (char*)malloc((pathLen+1)*sizeof(char));
        strncpy(path, argv[argc-3], pathLen+1);

        fileLen = strlen(argv[argc-2]);
        files[0] = (char*)malloc((fileLen+1)*sizeof(char));
        strncpy(files[0], argv[argc-2], fileLen+1);

        fileLen = strlen(argv[argc-1]);
        files[1] = (char*)malloc((fileLen+1)*sizeof(char));
        strncpy(files[1], argv[argc-1], fileLen+1);
        
        numberOfFiles = 2;
    } else {
        printf("Missing arguments!\n\n");
        printf("To compute distance of all fastq files from a directory use command:\n");
        printf("./gcBB <path_to_dir> --options\n\n");
        printf("To compute distance of two fastq files use command:\n");
        printf("./gcBB <path_to_dir> <file1> <file2> --options\n\n");
        exit(-1);
    }

    int systemTmp = system("mkdir tmp");
    if(systemTmp == -1){
        printf("Error creating tmp folder");
    }
    int systemResults = system("mkdir results");
    if(systemResults == -1){
        printf("Error creating results folder");
    }

    qsort(files, numberOfFiles, sizeof(char*), compareFiles);

    /******** Check PSUTIL ********/
    int result = system("python3 -c \"import psutil\" 2>/dev/null");
    if (result != 0) {
        printf("The 'psutil' library is NOT installed.\n");
        exit(-1);
    }

    /******** Compute external needed files ********/
    printf("=== PHASE 1 ===\n");
    printf("Start computing SA, BWT and LCP for all files\n");
    // Computes SA, BWT, LCP and DA from both files
    for(i = 0; i < numberOfFiles; i++){
        computeFile(path, files[i], memory);
    }

    printf("All needed arrays computed!\n");

    // Remove file format from the string
    for(i = 0; i < numberOfFiles; i++){
        char *ptr;
        ptr = strchr(files[i], '.');
        if (ptr != NULL)
            *ptr = '\0';
    }

    path = getPathDirName(path, pathLen);

    printf("Merging all pairs and computing document array (cda)\n");

    // Computes merge of files
    #if ALL_VS_ALL
        computeMergeFileAll(path, files, numberOfFiles, memory);
    #else
        for(i = 0; i < numberOfFiles; i++){
            for(j = i+1; j < numberOfFiles; j++){
                computeMergeFiles(path, files[i], files[j], memory);
            }
        }
    #endif

    printf("All arrays merged\n");

    // Similarity matrix based on expectation
    double **Dm = (double**)malloc(numberOfFiles*sizeof(double*));
    for(i = 0; i < numberOfFiles; i++)
        Dm[i] = (double*)malloc(numberOfFiles*sizeof(double));

    // Similarity matrix based on shannons entropy
    double **De = (double**)malloc(numberOfFiles*sizeof(double*));
    for(i = 0; i < numberOfFiles; i++)
        De[i] = (double*)malloc(numberOfFiles*sizeof(double));

    // Initialize matrixes
    for(i = 0; i < numberOfFiles; i++){
        for(j = 0; j < numberOfFiles; j++){
            Dm[i][j] = 0.0;
            De[i][j] = 0.0;
        }
    }

    printf("Start construction of colored BOSS and comparing genomes using BWSD for every pair\n");
    
    #if !ALL_VS_ALL
    for(i = 0; i < numberOfFiles; i++){
        for(j = i+1; j < numberOfFiles; j++){
            printf("=== PHASE 2 [%d,%d] ===\n", i, j);
    #else 
        printf("=== PHASE 2 ===\n");
    #endif
    char mergeBWTFile[FILE_PATH];
    char mergeLCPFile[FILE_PATH];
    char mergeDAFile[FILE_PATH];
    char mergeSLFile[FILE_PATH];

    #if !ALL_VS_ALL
        snprintf(mergeBWTFile, FILE_PATH, "tmp/merge.%s-%s.bwt", files[i], files[j]);
        snprintf(mergeLCPFile, FILE_PATH, "tmp/merge.%s-%s.2.lcp", files[i], files[j]);
        snprintf(mergeDAFile, FILE_PATH, "tmp/merge.%s-%s.1.cda", files[i], files[j]);
        snprintf(mergeSLFile, FILE_PATH, "tmp/merge.%s-%s.2.sl", files[i], files[j]);
    #else
        snprintf(mergeBWTFile, FILE_PATH, "tmp/merge.%s.bwt", path);
        snprintf(mergeLCPFile, FILE_PATH, "tmp/merge.%s.2.lcp", path);
        snprintf(mergeDAFile, FILE_PATH, "tmp/merge.%s.1.cda", path);
        snprintf(mergeSLFile, FILE_PATH, "tmp/merge.%s.2.sl", path);
    #endif

    FILE *mergeBWT = fopen(mergeBWTFile, "r");
    FILE *mergeLCP = fopen(mergeLCPFile, "rb");
    FILE *mergeDA = fopen(mergeDAFile, "rb");
    FILE *mergeSL = fopen(mergeSLFile, "rb");

    fseek(mergeBWT, 0, SEEK_END);
    size_t n = ftell(mergeBWT);
    rewind(mergeBWT);

    /******** Construct BOSS representation ********/
    #if !ALL_VS_ALL
        int samples = 2;
    #else
        int samples = numberOfFiles;
    #endif

    #if !ALL_VS_ALL
        bossConstruction(mergeLCP, mergeDA, mergeBWT, mergeSL, n, k, samples, memory, files[i], files[j], printBoss);
    #else
        bossConstruction(mergeLCP, mergeDA, mergeBWT, mergeSL, n, k, samples, memory, path, NULL, printBoss);
    #endif

    fclose(mergeBWT);
    fclose(mergeLCP);
    fclose(mergeDA);
    fclose(mergeSL);

    #if ALL_VS_ALL
        printf("=== PHASE 3 ===\n");
        #if COVERAGE
            bwsdAll(path, numberOfFiles, k, memory, Dm, De);
        #else 
            bwsdAll(path, numberOfFiles, k, memory, Dm, De);
        #endif
    #endif

    #if ALL_VS_ALL
    printf("For more details check file: results/%s_k_%d.info\n", path, k);
    #endif

    #if !ALL_VS_ALL
        printf("=== PHASE 3 [%d,%d] ===\n", i, j);
        double expectation, entropy;
        expectation = entropy = 0.0;
        bwsd(files[i], files[j], k, &expectation, &entropy, memory, printBoss, 0, 1);            
        Dm[j][i] = expectation;
        De[j][i] = entropy;
    #endif

    #if !ALL_VS_ALL
        printf("For more details check file: results/%s-%s_k_%d.info\n", files[i], files[j], k);
        }
    }
    #endif

    #if !ALL_VS_ALL
        printf("All genome pairs constructed and compared\n\n");
    #else
        printf("All genomes constructed and compared\n\n");
    #endif

    #if ALL_VS_ALL
    if(!printBoss){
        char colorFileName[FILE_PATH];
        char summarizedLCPFileName[FILE_PATH];
        char summarizedSLFileName[FILE_PATH];
        char coverageFileName[FILE_PATH];
        snprintf(colorFileName, FILE_PATH, "results/%s_k_%d.2.colors", path, k);
        snprintf(summarizedLCPFileName, FILE_PATH, "results/%s_k_%d.2.summarizedLCP", path, k);
        snprintf(summarizedSLFileName, FILE_PATH, "results/%s_k_%d.2.summarizedSL", path, k);
        snprintf(coverageFileName, FILE_PATH, "results/%s_k_%d.4.coverage", path, k);
        remove(colorFileName);
        remove(summarizedLCPFileName);
        remove(summarizedSLFileName);
        remove(coverageFileName);
    }
    #endif

    // Print BWSD results in files .dmat and .nhx
    printDistanceMatrixes(Dm, De, files, numberOfFiles, path, k);

    printf("All distance matrixes and newick files can be found in results folder\n");

    // Free variables
    for(i = 0; i < 512; i++) free(files[i]);
    free(files);

    for(i = 0; i < numberOfFiles; i++) free(Dm[i]);
    free(Dm);

    for(i = 0; i < numberOfFiles; i++) free(De[i]);
    free(De);

    free(path);
}

