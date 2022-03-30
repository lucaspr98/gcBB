#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <dirent.h>
#include <unistd.h>
#include <libgen.h>

#include "bwsd.h"
#include "boss.h"

#define FILE_PATH 1024

#ifndef COVERAGE
	#define COVERAGE 0
#endif

void compute_file(char *path, char *file);

void compute_merge_files(char *path, char *file1, char *file2);

void print_distance_matrixes(double **Dm, double **De, char **files, int files_n, char *path, int k, char coverage_type);

void compute_newick_files(char* dmat);

int compare_files(const void *element1, const void *element2) {
    const char **file_1 = (const char **)element1;
    const char **file_2 = (const char **)element2;
    return strcmp(*file_1, *file_2);
}

char drosophilas_files[11][64] = {
    "ananassae",
    "erecta",
    "melanogaster",
    "mojavensis",
    "persimilis",
    "pseudoobscura",
    "schellia",
    "simulans",
    "virilis",
    "willistoni",
    "yakuba"
};

int main(int argc, char *argv[]){
    int i, j;
    char **files = (char**)calloc(16, 64*sizeof(char*));
    int k = 30;
    int files_n = 11;
    char *path;
    int path_len;
    int opt;
    int memory = 1000;
    char coverage_type = 'a';
    int printBoss = 0;

    /******** Check arguments ********/

    int valid_opts = 0;
    while ((opt = getopt (argc, argv, "pk:m:c:")) != -1){
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
            case 'c':
                valid_opts += 2;
                coverage_type = *optarg;
                break;
            case '?':
                if(opt == 'k')
                    fprintf (stderr, "Option -%c requires a integer value.\n", opt);
                else if(opt == 'm')
                    fprintf (stderr, "Option -%c requires a integer value.\n", opt);
                else if(opt == 'c')
                    fprintf (stderr, "Option -%c requires one of the following options:\n'a': always apply coverage\n'e': apply coverage on same k-mers of distinct genomes\n'd': apply coverage in distinct k-mers", opt);
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
        for(i = 0; i < files_n; i++){
            files[i] = malloc((strlen(drosophilas_files)+1)*sizeof(char));
            strcpy(files[i], drosophilas_files[i]); 
        }	 
        DIR *folder;
        struct dirent *entry;
        int len;
        path_len = strlen(argv[argc-1]);
        path = (char*)malloc((path_len+1)*sizeof(char));
        strcpy(path, argv[argc-1]);

        // folder = opendir(argv[argc-1]);
        // if(folder == NULL){
        //     fprintf (stderr, "Unable to read directory %s\n", argv[argc-1]);
        //     exit(-1);
        // }
        // while((entry=readdir(folder)) != NULL){
        //     char *isFastq = strstr(entry->d_name, ".fastq");
        //     char *isFasta = strstr(entry->d_name, ".fasta");

        //     if((isFastq && strlen(isFastq) == 6) || (isFasta && strlen(isFasta) == 6)){
        //         len = strlen(entry->d_name)+1;
        //         files[files_n] = (char*)malloc((path_len+len+2)*sizeof(char));

        //         strcpy(files[files_n], entry->d_name);
        //         files_n++;
        //     }
        // }

        // closedir(folder);
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

    qsort(files, files_n, sizeof(char*), compare_files);

    /******** Compute external needed files ********/

    // Computes SA, BWT, LCP and DA from both files
    // for(i = 0; i < files_n; i++){
    //    compute_file(path, files[i]);
    // }

    // Remove file format from the string
    // for(i = 0; i < files_n; i++){
    //    char *ptr;
    //    ptr = strchr(files[i], '.');
    //    if (ptr != NULL)
    //        *ptr = '\0';
    // }

    // Computes merge of files
    // for(i = 0; i < files_n; i++){
    //    for(j = i+1; j < files_n; j++){
    //        compute_merge_files(path, files[i], files[j]);
    //    }
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
            Dm[i][j] = 0.0;
            De[i][j] = 0.0;
        }
    }

    for(i = 0; i < files_n-1; i++){
        for(j = i+1; j < files_n; j++){
		if(i == 1 || j == 1){
			printf("%s-%s\n", files[i], files[j]);
		
                
		char mergeBWTFile[FILE_PATH];
                char mergeLCPFile[FILE_PATH];
                char mergeDAFile[FILE_PATH];
                char mergeSLFile[FILE_PATH];

                sprintf(mergeBWTFile, "../drosophilas_merge/merge.%s-%s.bwt", files[i], files[j]);
                sprintf(mergeLCPFile, "../drosophilas_merge/merge.%s-%s.2.lcp", files[i], files[j]);
                sprintf(mergeDAFile, "../drosophilas_merge/merge.%s-%s.1.cda", files[i], files[j]);
                sprintf(mergeSLFile, "../drosophilas_merge/merge.%s-%s.2.sl", files[i], files[j]);

                FILE *mergeBWT = fopen(mergeBWTFile, "r");
                FILE *mergeLCP = fopen(mergeLCPFile, "rb");
                FILE *mergeDA = fopen(mergeDAFile, "rb");
                FILE *mergeSL = fopen(mergeSLFile, "rb");

                fseek(mergeBWT, 0, SEEK_END);
                size_t n = ftell(mergeBWT);
                rewind(mergeBWT);
  

              /******** Construct BOSS representation ********/
	        int samples = 2;

                size_t total_coverage = 0;

                size_t boss_len = boss_construction(mergeLCP, mergeDA, mergeBWT, mergeSL, n, k, samples, memory, files[i], files[j], printBoss, coverage_type, &total_coverage);

                fclose(mergeBWT);
                fclose(mergeLCP);
                fclose(mergeDA);

                /******** Compute BWSD ********/
                double expectation, entropy;
                expectation = entropy = 0.0;

                bwsd(files[i], files[j], boss_len, k, &expectation, &entropy, memory, printBoss, coverage_type, total_coverage);

                Dm[j][i] = expectation;
                De[j][i] = entropy;
    		}
  	  }
    }    

    // Print BWSD results in files .dmat and .nhx
    print_distance_matrixes(Dm, De, files, files_n, path, k, coverage_type);

    // Free variables
    for(i = 0; i < 32; i++) free(files[i]);
    free(files);

    for(i = 0; i < files_n; i++) free(Dm[i]);
    free(Dm);

    for(i = 0; i < files_n; i++) free(De[i]);
    free(De);

    free(path);


}

void compute_file(char *path, char *file){
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
        sprintf(eGap, "egap/eGap %s%s -m 12000 --em --rev --lcp  --sl --slbytes 2 -o tmp/%s", path, buff, file);
        system(eGap);    
    } else
        fclose(tmp);
}

void compute_merge_files(char *path, char *file1, char *file2){
    int len1 = strlen(file1); 
    int len2 = strlen(file2);

    char output[len1+len2+12];
    sprintf(output, "tmp/merge.%s-%s.bwt", file1, file2);
    FILE *tmp = fopen(output, "r");
    if(!tmp){
        char eGapMerge[FILE_PATH];
        sprintf(eGapMerge, "egap/eGap -m 12000 --em --bwt --lcp --cda --cbytes 1 --sl --slbytes 2 --rev tmp/%s.bwt tmp/%s.bwt -o tmp/merge.%s-%s", file1, file2, file1, file2);
        system(eGapMerge);    
    } else 
        fclose(tmp);
}

void print_distance_matrixes(double **Dm, double **De, char **files, int files_n, char *path, int k, char coverage_type){
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
            char coverage_arg[FILE_PATH];
            sprintf(coverage_arg, "_coverage_%c", coverage_type);
            strcat(expectationDmat, coverage_arg);
            strcat(entropyDmat, coverage_arg);
        #endif

        char extension[FILE_PATH];
        sprintf(extension, "_k_%d.dmat", k);
        strcat(expectationDmat, extension);
        strcat(entropyDmat, extension);

    } else {
        sprintf(expectationDmat, "results/%s-%s_expectation", files[0], files[1]);            
        sprintf(entropyDmat, "results/%s-%s_entropy", files[0], files[1]);            

        #if COVERAGE
            char coverage_arg[FILE_PATH];
            sprintf(coverage_arg, "_coverage_%c", coverage_type);
            strcat(expectationDmat, coverage_arg);
            strcat(entropyDmat, coverage_arg);
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
    
    compute_newick_files(expectationDmat);
    compute_newick_files(entropyDmat);
}


void compute_newick_files(char *dmat){
    char dmat_newick[FILE_PATH];

    char *ptr = strchr(dmat, '.');
    if (ptr != NULL)
        *ptr = '\0';

    sprintf(dmat_newick, "utils/nj -i %s.dmat -n %s.nhx", dmat, dmat);

    system(dmat_newick);
}

