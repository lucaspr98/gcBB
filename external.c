#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <dirent.h>
#include <unistd.h>
#include <libgen.h>
#include "external.h"

#define FILE_PATH 1024

void computeNewickFiles(char *dmat){
    char dmat_newick[FILE_PATH];

    char *ptr = strchr(dmat, '.');
    if (ptr != NULL)
        *ptr = '\0';

    sprintf(dmat_newick, "utils/nj -i %s.dmat -n %s.nhx", dmat, dmat);

    system(dmat_newick);
}

void computeFile(char *path, char *file, int memory){
    int len = strlen(file);
    char buff[len+1];
    strcpy(buff, file); 

    char *ptr = strchr(file, '.');
    if(ptr != NULL)
        *ptr = '\0';

    char output[len+10];
    sprintf(output, "tmp/%s.bwt", file);
    FILE *tmp = fopen(output, "r");
    if(!tmp){
        char eGap[FILE_PATH];
        sprintf(eGap, "egap/eGap %s%s -m %d --em --rev --lcp  --sl --slbytes 2 -o tmp/%s", path, buff, memory, file);
        system(eGap);    
    } else {
        printf("%s files already computed!\n", file);
        fclose(tmp);
    }
}

void computeMergeFileAll(char *path, char **files, int numberOfFiles, int memory){
    char output[strlen(path)+12];
    sprintf(output, "tmp/merge.%s.bwt", path);
    FILE *tmp = fopen(output, "r");
    if(!tmp){
        char eGapMerge[FILE_PATH];
        sprintf(eGapMerge, "egap/eGap -m %d --em --bwt --lcp --cda --cbytes 1 --sl --slbytes 2 ", memory);
        for(int i = 0; i < numberOfFiles; i++)
            sprintf(eGapMerge + strlen(eGapMerge), " tmp/%s.bwt ", files[i]);
        sprintf(eGapMerge + strlen(eGapMerge), " -o tmp/merge.%s", path);
        printf("%s\n", eGapMerge);
        system(eGapMerge);    
    } else {
        printf("%s merge file already computed!\n", path);
        fclose(tmp);
    }
}

void computeMergeFiles(char *path, char *file1, char *file2, int memory){
    int len1 = strlen(file1); 
    int len2 = strlen(file2);

    char output[len1+len2+12];
    sprintf(output, "tmp/merge.%s-%s.bwt", file1, file2);
    FILE *tmp = fopen(output, "r");
    if(!tmp){
        char eGapMerge[FILE_PATH];
        sprintf(eGapMerge, "egap/eGap -m %d --em --bwt --lcp --cda --cbytes 1 --sl --slbytes 2 --rev tmp/%s.bwt tmp/%s.bwt -o tmp/merge.%s-%s", memory, file1, file2, file1, file2);
        system(eGapMerge);    
    } else {
        printf("%s-%s merge files already computed!\n", file1, file2);
        fclose(tmp);
    }
}

void printDistanceMatrixes(double **Dm, double **De, char **files, int files_n, char *path, int k){
    int i,j;
    char *ptr;

    int len = strlen(path);
    char folder[len];
    strcpy(folder, basename(path));

    char expectationDmat[FILE_PATH];
    char entropyDmat[FILE_PATH];

    if(files_n > 2){
        ptr = strchr(path, '/');
        if (ptr != NULL)
            *ptr = '\0';

        sprintf(expectationDmat, "results/%s_expectation", folder);
        sprintf(entropyDmat, "results/%s_entropy", folder);

        #if COVERAGE
            strcat(expectationDmat, "_coverage");
            strcat(entropyDmat, "_coverage");
        #endif
        #if ALL_VS_ALL
            strcat(expectationDmat, "_all");
            strcat(entropyDmat, "_all");
        #endif

        char extension[FILE_PATH];
        sprintf(extension, "_k_%d.dmat", k);
        strcat(expectationDmat, extension);
        strcat(entropyDmat, extension);

    } else {
        sprintf(expectationDmat, "results/%s_expectation", path);            
        sprintf(entropyDmat, "results/%s_entropy",  path);            

        #if COVERAGE
            strcat(expectationDmat, "_coverage");
            strcat(entropyDmat, "_coverage");
        #endif
        #if ALL_VS_ALL
            strcat(expectationDmat, "_all");
            strcat(entropyDmat, "_all");
        #endif

        char extension[FILE_PATH];
        sprintf(extension, "_k_%d.dmat", k);
        strcat(expectationDmat, extension);
        strcat(entropyDmat, extension);
    }

    FILE *expectationDmatFile = fopen(expectationDmat, "w");
    FILE *entropyDmatFile = fopen(entropyDmat, "w");

    fprintf(expectationDmatFile, "[size]\n%d\n", files_n);

    fprintf(expectationDmatFile, "[labels]\n");
    for(i = 0; i < files_n; i++){
        fprintf(expectationDmatFile, "%s ", files[i]);
    }

    fprintf(expectationDmatFile, "\n");

    fprintf(expectationDmatFile, "[distances]\n");
    for(i = 1; i < files_n; i++){
        for(j = 0; j < i; j++){
            fprintf(expectationDmatFile, "%lf\t", Dm[i][j]);
        }
        fprintf(expectationDmatFile, "\n");
    }

    fprintf(entropyDmatFile, "[size]\n%d\n", files_n);

    fprintf(entropyDmatFile, "[labels]\n");
    for(i = 0; i < files_n; i++){
        fprintf(entropyDmatFile, "%s ", files[i]);
    }

    fprintf(entropyDmatFile, "\n");

    fprintf(entropyDmatFile, "[distances]\n");
    for(i = 1; i < files_n; i++){
        for(j = 0; j < i; j++){
            fprintf(entropyDmatFile, "%lf\t", De[i][j]);
        }
        fprintf(entropyDmatFile, "\n");
    }


    fclose(expectationDmatFile);
    fclose(entropyDmatFile);
    
    printf("entropy newick file construction:\n");
    computeNewickFiles(entropyDmat);
    printf("expectation newick file construction:\n");
    computeNewickFiles(expectationDmat);
}

FILE* getInfoFile(char* file1, char* file2, int k, int update){
    char info[FILE_PATH];

    #if ALL_VS_ALL
    sprintf(info, "results/%s_k_%d", file1, k);
    #else
    sprintf(info, "results/%s-%s_k_%d", file1, file2, k);   
    #endif

    #if COVERAGE
    strcat(info, "_cov");
    #endif

    #if ALL_VS_ALL
    strcat(info, "_all");
    #endif

    #if DEBUG
    strcat(info, ".info.debug");
    #else
    strcat(info, ".info");
    #endif

    if(update){
        return fopen(info, "a+");
    } else
        return fopen(info, "w");
}