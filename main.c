#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <dirent.h>
#include <unistd.h>

#include "bwsd.h"
#include "boss.h"

void compute_file(char *path, char *file, int k);

void compute_merge_files(char *path, char *file1, char *file2, int k);

int main(int argc, char *argv[]){
    int i, j, m;
    char **files = (char**)malloc(32*sizeof(char*));
    int k = 30;
    int files_n = 0;
    char *path;
    int path_len;
    int opt;
    int memory = 1000;

    /******** Check arguments ********/

    int valid_opts = 0;
    while ((opt = getopt (argc, argv, "k:m:")) != -1){
        switch (opt){
            case 'k':
                valid_opts++;
                k = atoi(optarg);
                break;
            case 'm':
                valid_opts++;
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

    if(argc-2*valid_opts == 2){
        DIR *folder;
        struct dirent *entry;
        int len;
        path_len = strlen(argv[argc-1]);
        path = (char*)malloc(path_len*sizeof(char));
        strcpy(path, argv[argc-1]);

        folder = opendir(argv[argc-1]);
        if(folder == NULL){
            fprintf (stderr, "Unable to read directory %s\n", argv[argc-1]);
            exit(-1);
        }
        while((entry=readdir(folder)) != NULL){
            char *isFastq = strstr(entry->d_name, ".fastq");
            char *isFasta = strstr(entry->d_name, ".fasta");

            if((isFastq && strlen(isFastq) == 6) || (isFasta && strlen(isFasta) == 6)){
                printf("%s\n", entry->d_name);
                len = strlen(entry->d_name)+1;
                files[files_n] = (char*)malloc((path_len+len+2)*sizeof(char));

                strcpy(files[files_n], entry->d_name);
                files_n++;
            }
        }

        closedir(folder);
    } else if(argc-2*valid_opts == 4){
        int file_len;
        path_len = strlen(argv[argc-3]);
        path = (char*)malloc(path_len*sizeof(char));
        strcpy(path, argv[argc-3]);

        file_len = strlen(argv[argc-2]);
        files[0] = (char*)malloc((file_len+1)*sizeof(char));
        strcpy(files[0], argv[argc-2]);

        file_len = strlen(argv[argc-1]);
        files[1] = (char*)malloc((file_len+1)*sizeof(char));
        strcpy(files[1], argv[argc-1]);
        
        files_n = 2;
    } else {
        printf("Missing arguments!\n\n");
        printf("To compute distance of all fastq files from a directory use command:\n");
        printf("./gcBB <path_to_dir> --options\n\n");
        printf("To compute distance of two fastq files use command:\n");
        printf("./gcBB <path_to_dir> <file1> <file2> --options\n\n");
        exit(-1);
    }

    system("mkdir tmp");
    system("mkdir results");

    /******** Compute external needed files ********/

    // Computes SA, BWT, LCP and DA from both files
    // for(i = 0; i < files_n; i++){
    //     compute_file(path, files[i], k);
    // }

    // Remove file format from the string
    for(i = 0; i < files_n; i++){
        int len = strlen(files[i]);
        char *buff = malloc(len*sizeof(char));
        strncpy(buff, files[i], len-6);
        strcpy(files[i], buff);
        printf("%s\n", files[i]);
    }
    
    // Computes merge of files
    // for(i = 0; i < files_n; i++){
    //     for(j = i+1; j < files_n; j++){
    //         compute_merge_files(path, files[i], files[j], k);
    //     }
    // }

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
                char mergeBWTFile[128];
                char mergeLCPFile[128];
                char mergeDAFile[128];

                sprintf(mergeBWTFile, "tmp/merge.%s-%s.bwt", files[i], files[j]);
                sprintf(mergeLCPFile, "tmp/merge.%s-%s.2.lcp", files[i], files[j]);
                sprintf(mergeDAFile, "tmp/merge.%s-%s.4.da", files[i], files[j]);

                FILE *mergeBWT = fopen(mergeBWTFile, "r");
                FILE *mergeLCP = fopen(mergeLCPFile, "rb");
                FILE *mergeDA = fopen(mergeDAFile, "rb");

                fseek(mergeBWT, 0, SEEK_END);
                size_t n = ftell(mergeBWT);
                rewind(mergeBWT);

                /******** Construct BOSS representation ********/

                // BOSS construction needed variables
                // char *BWT;
                // short *LCP;
                // int *DA;
                int *last, *Wm, *colors;
                char *W;
                int C[255] = {0};
                int samples = 2;

                // Coverage information variables
                int total_coverage = 0;
                short *reduced_LCP;
                int *coverage;

                reduced_LCP = (short*)malloc(n*sizeof(short));
                coverage = (int*)malloc(n*sizeof(int));

                // Initialize variables needed to construct BOSS
                // BWT = (char*)malloc(n*sizeof(char));
                // LCP = (short*)malloc(n*sizeof(short));
                // DA = (int*)malloc(n*sizeof(int));

                // fread(BWT, sizeof(char), n, mergeBWT);
                // fread(LCP, sizeof(short), n, mergeLCP);
                // fread(DA, sizeof(int), n, mergeDA);
                // for(m = 0; m < n; m++){
                //     if(BWT[m] == 0)
                //         BWT[m] = '$';
                //     else
                //         BWT[m] = toupper(BWT[m]);

                //     coverage[m] = 1;
                // }

                char docsFile[128];
                sprintf(docsFile, "tmp/%s.docs", files[i]);
                FILE *docs = fopen(docsFile, "r");
                size_t docsSeparator;
                fread(&docsSeparator, 8, 1, docs);

                // for(m = 0; m < n; m++){
                //     if(DA[m] < docsSeparator)
                //         DA[m] = 0;
                //     else
                //         DA[m] = 1;
                // }

                // Initialize BOSS variables
                last = (int*)malloc(n*sizeof(int));
                Wm = (int*)malloc(n*sizeof(int));
                colors = (int*)malloc(n*sizeof(int));
                W = (char*)malloc(n*sizeof(char));

                int boss_len = boss_construction(mergeLCP, mergeDA, mergeBWT, C, last, W, Wm, colors, n, k, samples, reduced_LCP, coverage, &total_coverage, docsSeparator, memory);

                // Print BOSS result
                // print_boss_result(boss_len, files[i], files[j], C, last, W, Wm, colors, reduced_LCP, coverage, total_coverage);
                
                // free(LCP);free(BWT);
                free(last);free(W);free(Wm);free(colors);


                fclose(mergeBWT);
                fclose(mergeLCP);

                /******** Compute BWSD ********/
                double expectation, entropy;
                expectation = entropy = 0;

                // bwsd(mergeDA, reduced_LCP, coverage, boss_len, k, &expectation, &entropy, docsSeparator, memory);


                fclose(mergeDA);
        
                Dm[i][j] = expectation;
                De[i][j] = entropy;

                // free(DA); 
                free(reduced_LCP);
            } else if (j == i){
                Dm[i][j] = 0;
                De[i][j] = 0;
            }
        }
    }    

    // Print BWSD results
    // print_bwsd_matrixes(Dm, De, files, files_n);

    // system("rm -rf tmp");
}

void compute_file(char *path, char *file, int k){
    char eGap[256];
    int len = strlen(file);
    char *buff = malloc(len*sizeof(char));
    strncpy(buff, file, len-6);

    sprintf(eGap, "egap/eGap %s%s -m 8192 --em --rev --da  --lcp -o tmp/%s", path, file, buff);

    system(eGap);
}

void compute_merge_files(char *path, char *file1, char *file2, int k){
    char eGapMerge[256];

    sprintf(eGapMerge, "egap/eGap -m 8192 --em --bwt --trlcp %d --da --rev tmp/%s.bwt tmp/%s.bwt -o tmp/merge.%s-%s", k+1, file1, file2, file1, file2);
    
    system(eGapMerge);
}