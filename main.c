#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <dirent.h>
#include <unistd.h>
#include <libgen.h>

#include "bwsd.h"
#include "boss.h"

void compute_file(char *path, char *file, int k);

void compute_merge_files(char *path, char *file1, char *file2, int k);

int main(int argc, char *argv[]){
    int i, j, m;
    char **files = (char**)calloc(10, 32*sizeof(char*));
    int k = 30;
    int files_n = 0;
    char *path;
    int path_len;
    int opt;
    int memory = 1000;
    int printBoss = 0;

    /******** Check arguments ********/

    int valid_opts = 0;
    while ((opt = getopt (argc, argv, "pk:m:")) != -1){
        switch (opt){
            case 'p':
                valid_opts+=1;
                printBoss = 1;
                break;
            case 'k':
                valid_opts += 2;
                k = atoi(optarg);
                break;
            case 'm':
                valid_opts += 2;
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

    if(argc-valid_opts == 2){
        DIR *folder;
        struct dirent *entry;
        int len;
        path_len = strlen(argv[argc-1]);
        path = (char*)malloc((path_len+1)*sizeof(char));
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
    } else if(argc-valid_opts == 4){
        int file_len;
        path_len = strlen(argv[argc-3]);
        path = (char*)malloc((path_len+1)*sizeof(char));
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
    for(i = 0; i < files_n; i++){
        compute_file(path, files[i], k);
    }

    // Remove file format from the string
    for(i = 0; i < files_n; i++){
        int len = strlen(files[i]);
        char *ptr;
        ptr = strchr(files[i], '.');
        if (ptr != NULL)
            *ptr = '\0';
    }
    
    // Computes merge of files
    for(i = 0; i < files_n; i++){
        for(j = i+1; j < files_n; j++){
            compute_merge_files(path, files[i], files[j], k);
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
                char mergeBWTFile[128];
                char mergeLCPFile[128];
                char mergeDAFile[128];

                sprintf(mergeBWTFile, "tmp/merge.%s-%s.bwt", files[i], files[j]);
                sprintf(mergeLCPFile, "tmp/merge.%s-%s.2.lcp", files[i], files[j]);
                sprintf(mergeDAFile, "tmp/merge.%s-%s.1.cda", files[i], files[j]);

                FILE *mergeBWT = fopen(mergeBWTFile, "r");
                FILE *mergeLCP = fopen(mergeLCPFile, "rb");
                FILE *mergeDA = fopen(mergeDAFile, "rb");

                fseek(mergeBWT, 0, SEEK_END);
                size_t n = ftell(mergeBWT);
                rewind(mergeBWT);

                /******** Construct BOSS representation ********/

                int samples = 2;

                size_t boss_len = boss_construction(mergeLCP, mergeDA, mergeBWT, n, k, samples, memory, files[i], files[j]);

                fclose(mergeBWT);
                fclose(mergeLCP);
                fclose(mergeDA);

                /******** Compute BWSD ********/
                double expectation, entropy;
                expectation = entropy = 0;

                bwsd(files[i], files[j], boss_len, k, &expectation, &entropy, memory);

                Dm[i][j] = expectation;
                De[i][j] = entropy;
            } else if (j == i){
                Dm[i][j] = 0;
                De[i][j] = 0;
            }
        }
    }    

    // Print BWSD results
    print_bwsd_matrixes(Dm, De, files, files_n, path);

    // Free variables
    for(i = 0; i < 32; i++) free(files[i]);
    free(files);

    for(i = 0; i < files_n; i++) free(Dm[i]);
    free(Dm);

    for(i = 0; i < files_n; i++) free(De[i]);
    free(De);

    free(path);

    // system("rm -rf tmp");
}

void compute_file(char *path, char *file, int k){
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
        char eGap[256];
        sprintf(eGap, "egap/eGap %s%s -m 12288 --em --rev --lcp -o tmp/%s", path, buff, file);
        system(eGap);    
    } else
        fclose(tmp);
}

void compute_merge_files(char *path, char *file1, char *file2, int k){
    int len1 = strlen(file1); 
    int len2 = strlen(file2);
    int tmpSize = len1+len2+4;

    char output[len1+len2+12];
    sprintf(output, "tmp/merge.%s-%s.bwt", file1, file2);
    FILE *tmp = fopen(output, "r");
    if(!tmp){
        char eGapMerge[256];
        sprintf(eGapMerge, "egap/eGap -m 12288 --em --bwt --trlcp %d --cda --cbytes 1 --rev tmp/%s.bwt tmp/%s.bwt -o tmp/merge.%s-%s", k, file1, file2, file1, file2);
        system(eGapMerge);    
    } else 
        fclose(tmp);
}