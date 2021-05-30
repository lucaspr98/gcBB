#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <dirent.h>

#include "bwsd.h"
#include "boss.h"

void compute_file(char *path, char *file, int k);

void compute_merge_files(char *path, char *file1, char *file2, int id1, int id2, int k);

int main(int argc, char *argv[]){
    int i, j, m;
    char **files = (char**)malloc(20*sizeof(char*));
    int k;
    int files_n = 0;
    char *path;
    int path_len;

    /******** Check arguments ********/

    if (argc == 3){
        DIR *folder;
        struct dirent *entry;
        int len;
        path_len = strlen(argv[1]);
        path = (char*)malloc(path_len*sizeof(char));
        strcpy(path, argv[1]);

        folder = opendir(argv[1]);
        if(folder == NULL){
            perror("Unable to read directory");
        }
        while((entry=readdir(folder)) != NULL){
            char *isfastq = strstr(entry->d_name, ".fastq");

            if(isfastq && strlen(isfastq) == 6){
                printf("%s\n", entry->d_name);
                len = strlen(entry->d_name)+1;
                files[files_n] = (char*)malloc((path_len+len+2)*sizeof(char));

                strcpy(files[files_n], entry->d_name);
                files_n++;
            }
        }

        closedir(folder);

        k = atoi(argv[2]);
    } else if (argc == 5){
        int file_len;
        path_len = strlen(argv[1]);
        path = (char*)malloc(path_len*sizeof(char));
        strcpy(path, argv[1]);

        file_len = strlen(argv[2]);
        files[0] = (char*)malloc((file_len+1)*sizeof(char));
        file_len = strlen(argv[3]);
        files[1] = (char*)malloc((file_len+1)*sizeof(char));

        strcpy(files[0], argv[2]);
        strcpy(files[1], argv[3]);
        
        files_n = 2;

        k = atoi(argv[4]);
    } else {
        printf("Missing arguments!\n\n");
        printf("To compute distance of all fastq files from a directory use command:\n");
        printf("./gcBB <path_to_dir> <k>\n\n");
        printf("To compute distance of two fastq files use command:\n");
        printf("./gcBB <path_to_dir> <file1.fastq> <file2.fastq> <k>\n\n");
        exit(-1);
    }    

    system("mkdir tmp");
    system("mkdir results");

    /******** Compute external needed files ********/

    // Computes SA, BWT, LCP and DA from both files
    for(i = 0; i < files_n; i++){
        compute_file(path, files[i], k);
    }
    
    // Computes merge of files
    for(i = 0; i < files_n; i++){
        for(j = i+1; j < files_n; j++){
            compute_merge_files(path, files[i], files[j], i+1, j+1, k);
        }
    }

    // Similarity matrix based on expectation
    double **Dm = (double**)malloc(files_n*sizeof(double*));
    for(i = 0; i < files_n; i++)
        Dm[i] = (double*)malloc(files_n*sizeof(double));

    // Similarity matrix based on shannons entropy
    double **De = (double**)malloc(files_n*sizeof(double*));
    for(i = 0; i < files_n; i++)
        De[i] = (double*)malloc(files_n*sizeof(double));

    // Initialize matrixes
    for(i = 0; i < files_n; i++){
        for(j = 0; j < files_n; j++){
            Dm[i][j] = 0;
            De[i][j] = 0;
        }
    }

    for(i = 0; i < files_n; i++){
        for(j = 0; j < files_n; j++){
            if(j > i){
                char mergeBWTFile[32];
                char mergeLCPFile[32];
                char mergeDAFile[32];

                sprintf(mergeBWTFile, "tmp/merge.%d-%d.bwt", i+1, j+1);
                sprintf(mergeLCPFile, "tmp/merge.%d-%d.4.lcp", i+1, j+1);
                sprintf(mergeDAFile, "tmp/merge.%d-%d.4.da", i+1, j+1);

                FILE *mergeBWT = fopen(mergeBWTFile, "r");
                FILE *mergeLCP = fopen(mergeLCPFile, "rb");
                FILE *mergeDA = fopen(mergeDAFile, "rb");

                fseek(mergeBWT, 0, SEEK_END);
                int n = ftell(mergeBWT);
                rewind(mergeBWT);

                /******** Construct BOSS representation ********/

                // BOSS construction needed variables
                char *BWT;
                int *LCP, *DA;
                int *last, *Wm, *colors;
                char *W;
                int C[255] = {0};
                int samples = 2;

                // Initialize variables needed to construct BOSS
                BWT = (char*)malloc(n*sizeof(char));
                LCP = (int*)malloc(n*sizeof(int));
                DA = (int*)malloc(n*sizeof(int));

                fread(BWT, sizeof(char), n, mergeBWT);
                fread(LCP, sizeof(int), n, mergeLCP);
                fread(DA, sizeof(int), n, mergeDA);
                for(m = 0; m < n; m++){
                    if(BWT[m] == 0)
                        BWT[m] = '$';
                    else
                        BWT[m] = toupper(BWT[m]);
                }

                fclose(mergeBWT);
                fclose(mergeLCP);
                fclose(mergeDA);

                // Initialize BOSS variables
                // memset(C, 0, sizeof(int)*255);
                last = (int*)malloc(n*sizeof(int));
                Wm = (int*)malloc(n*sizeof(int));
                colors = (int*)malloc(n*sizeof(int));
                W = (char*)malloc(n*sizeof(char));

                int boss_len = boss_construction(LCP, DA, BWT, C, last, W, Wm, colors, n, k, samples);

                // Print BOSS result
                print_boss_result(boss_len, i+1, j+1, files[i], files[j], C, last, W, Wm, colors);
                
                free(LCP);free(BWT);free(last);free(W);free(Wm);free(colors);

                /******** Compute BWSD ********/
                double expectation, entropy;
                expectation = entropy = 0;

                bwsd(DA, boss_len, &expectation, &entropy);
        
                Dm[i][j] = expectation;
                De[i][j] = entropy;

                free(DA);
            } else if (j == i){
                Dm[i][j] = 1;
                De[i][j] = 1;
            }
        }
    }    

    // Print BWSD results
    print_bwsd_matrixes(Dm, De, files, files_n);

    system("rm -rf tmp");
}

void compute_file(char *path, char *file, int k){
    char eGap[256];

    sprintf(eGap, "egap/eGap -m 4096 %s%s -o tmp/%s --trlcp %d --rev --lbytes 4 --da", path, file, file, k);
    
    system(eGap);
}

void compute_merge_files(char *path, char *file1, char *file2, int id1, int id2, int k){
    char eGapMerge[256];

    sprintf(eGapMerge, "egap/eGap -m 4096 --trlcp %d --rev --bwt --lbytes 4 --da -o tmp/merge.%d-%d tmp/%s.bwt tmp/%s.bwt", k, id1, id2, file1, file2);
    
    system(eGapMerge);
}