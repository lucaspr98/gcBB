#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <dirent.h>

void compute_file(char *path, char *file, int k);

void compute_merge_files(char *path, char *file1, char *file2, int id1, int id2, int k);

int boss_construction(int *LCP, int *DA, char *BWT, int *C, int *last, char *W, int *Wm, int *colors, int n, int k, int samples);

void Wi_sort(char *Wi, int *Wm,  int *colors, int start, int end);

void bwsd(int *DA, int n, double *expectation, double *entropy);

double bwsd_expectation(int *t, int s, int n);

double bwsd_shannon_entropy(int *t, int s, int n);

void print_boss_result(int boss_len, int id1, int id2, char *file1, char *file2, int *C, int *last, char *W, int *Wm, int *colors);

void print_bwsd_matrixes(double **Dm, double **De, char **files, int files_n);

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
        while((entry=readdir(folder))){
            char *isfastq = strstr(entry->d_name, ".fastq");
            if(isfastq && strlen(isfastq) == 6){
                len = strlen(entry->d_name)+1;
                files[files_n] = (char*)malloc((path_len+len+2)*sizeof(char));

                // char id = (files_n+1);
                // files[files_n][0] = '0'+id;
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

void Wi_sort(char *Wi, int *Wm, int *colors, int start, int end){
    int i;
    int range = end-start;
    char Wi_tmp[range];
    int Wm_tmp[range];
    int colors_tmp[range];
    int Wi_aux[255];
    int Wm_aux[255];
    int repetitive[255];

    memset(Wi_tmp, 0, sizeof(char)*(range));    
    memset(Wm_tmp, 0, sizeof(int)*(range)); 
    memset(colors_tmp, 0, sizeof(int)*(range)); 
    memset(Wi_aux, 0, sizeof(int)*255);
    memset(repetitive, 0, sizeof(int)*255);
    memset(Wm_aux, 0, sizeof(int)*255);

    for(i = start; i < end; i++){
        Wi_aux[Wi[i]]++;
        Wm_aux[Wi[i]] += Wm[i];
    }

    for(i = start; i < end; i++){
        if(Wi_aux[Wi[i]] > 1)
            repetitive[Wi[i]] = 1;
    }

    for(i = 1; i < 255; i++){
        Wi_aux[i] += Wi_aux[i-1];
    }

    for(i = start; i < end; i++){
        Wi_tmp[Wi_aux[Wi[i]]-1] = Wi[i];
        
        if(repetitive[Wi[i]] == 1 && Wm_aux[Wi[i]] == 0){
            Wm_tmp[Wi_aux[Wi[i]]-1] = 1;
        } else if (repetitive[Wi[i]] == 1 && Wm_aux[Wi[i]] > 0){
            Wm_tmp[Wi_aux[Wi[i]]-1] = 0;
            Wm_aux[Wi[i]]--;
        } else {
            Wm_tmp[Wi_aux[Wi[i]]-1] = Wm[i];
        }
        colors_tmp[Wi_aux[Wi[i]]-1] = colors[i];
        Wi_aux[Wi[i]]--;
    }

    for(i = start; i < end; i++){
        Wi[i] = Wi_tmp[i-start];
        Wm[i] = Wm_tmp[i-start];
        colors[i] = colors_tmp[i-start];
    }
}

int boss_construction(int *LCP, int *DA, char *BWT, int *C, int *last, char *W, int *Wm, int *colors, int n, int k, int samples){
    int i = 0; // iterates through Wi
    int j = 0; // auxiliary iterator  
    int bi = 0; // iterates through BWT and LCP
    int Wi_size = 0; 
    int W_freq[255]; // frequency of outgoing edges in a (k-1)-mer suffix range (detect W- = 1)
    memset(W_freq, 0, sizeof(int)*255);
    int Wi_freq[255]; // frequency of outgoing edges in a k-mer suffix range (detect same outgoing edge in a vertex)
    memset(Wi_freq, 0,  sizeof(int)*255);
    int DA_freq[samples][255]; // frequency of outgoing edges in a k-mer from a string collection (used to include same outgoing edge from distinct collections in BOSS representation)
    memset(DA_freq, 0,  sizeof(int)*samples*255);

    while(bi < n){
        // more than one outgoing edge of vertex i
        if(LCP[bi+1] >= k && bi != n-1){
            // since there is more than one outgoing edge, we don't need to check if BWT = $ or there is already BWT[bi] in Wi range
            if(BWT[i] != '$'){
                if(Wi_freq[BWT[bi]] == 0){
                    // Add values to BOSS representation
                    W[i] = BWT[bi];
                    colors[i] = DA[bi];
                    if(Wi_size == 0){
                        last[i] = 1;
                    } else {
                        last[i-1] = 0;
                        last[i] = 1;
                    }
                    if(W_freq[BWT[bi]] == 0){
                        Wm[i] = 1;
                    }
                    // Increment variables
                    C[BWT[bi]]++; W_freq[BWT[bi]]++; Wi_freq[BWT[bi]]++; DA_freq[DA[bi]][BWT[bi]]++; Wi_size++; i++; 
                } else {
                    // check if there is already outgoing edge labeled with BWT[bi] from DA[bi] leaving vertex i
                    if(DA_freq[DA[bi]][BWT[bi]] == 0){
                        W[i] = BWT[bi];
                        colors[i] = DA[bi];
                        if(Wi_size == 0){
                            last[i] = 1;
                        } else {
                            last[i-1] = 0;
                            last[i] = 1;
                        }
                        if(W_freq[BWT[bi]] == 0){
                            Wm[i] = 1;
                        }
                        C[BWT[bi]]++; W_freq[BWT[bi]]++; Wi_freq[BWT[bi]]++; DA_freq[DA[bi]][BWT[bi]]++; Wi_size++; i++; 
                    }
                }
            }
        } else {
            // just one outgoing edge of vertex i
            if(Wi_size == 0){
                W[i] = BWT[bi];
                colors[i] = DA[bi];
                last[i] = 1;
                if(W_freq[BWT[bi]] == 0){
                    Wm[i] = 1;
                }
                C[BWT[bi]]++; W_freq[BWT[bi]]++; i++;
            } 
            // last outgoing edge of vertex i
            else {
                // check if there is already outgoing edge labeled with BWT[bi] leaving vertex i
                if(Wi_freq[BWT[bi]] == 0){
                    W[i] = BWT[bi];
                    last[i-1] = 0;
                    last[i] = 1;
                    colors[i] = DA[bi];
                    if(W_freq[BWT[bi]] == 0){
                        Wm[i] = 1;
                    }
                    C[BWT[bi]]++; W_freq[BWT[bi]]++; Wi_size++; i++;
                } else {
                    // check if there is already outgoing edge labeled with BWT[bi] from DA[bi] leaving vertex i
                    if(DA_freq[DA[bi]][BWT[bi]] == 0){
                        W[i] = BWT[bi];
                        colors[i] = DA[bi];
                        if(Wi_size == 0){
                            last[i] = 1;
                        } else {
                            last[i-1] = 0;
                            last[i] = 1;
                        }
                        if(W_freq[BWT[bi]] == 0){
                            Wm[i] = 1;
                        }
                        C[BWT[bi]]++; W_freq[BWT[bi]]++; Wi_freq[BWT[bi]]++; DA_freq[DA[bi]][BWT[bi]]++; Wi_size++; i++; 
                    }
                }
                // sort outgoing edges of vertex i in lexigraphic order
                if(Wi_size > 1){
                    Wi_sort(W, Wm, colors, i-Wi_size, i);
                }
                // clean frequency variables of outgoing edges in Wi 
                memset(Wi_freq, 0, sizeof(int)*255);   
                memset(DA_freq, 0, sizeof(int)*samples*255);
            }
            // if next LCP value is smaller than k-1 we have a new (k-1)-mer to keep track, so we clean W_freq values
            if(LCP[bi+1] < k-1){
                memset(W_freq, 0,  sizeof(int)*255);
            }
            Wi_size = 0; 
        }
        bi++;
    }

    // fix C values
    C[1] = C['$'];
    C[2] = C['A'] + C[1];
    C[3] = C['C'] + C[2];
    C[4] = C['G'] + C[3];
    C[5] = C['N'] + C[4];
    C[0] = 0;

    return i;
};

void print_boss_result(int boss_len, int id1, int id2, char *file1, char *file2, int *C, int *last, char *W, int *Wm, int *colors){
    int i;
    char alphabet[6] = {'$', 'A', 'C', 'G', 'N', 'T'};
    char boss_result[64];
    
    sprintf(boss_result, "results/%d-%d.boss", id1, id2);
    
    FILE *boss_file = fopen(boss_result, "w");
                
    fprintf(boss_file, "Boss construction of %s and %s genomes merge:\n", file1, file2);
    fprintf(boss_file, "C array:\n");
    for(i = 0; i < 6; i++)
        fprintf(boss_file, "%c %d\n", alphabet[i], C[i]);
    fprintf(boss_file, "\n");
    
    fprintf(boss_file, "BOSS:\nlast\tW\tW-\tcolor\n");
    for(i = 0; i < boss_len; i++)
        fprintf(boss_file, "%d\t\t%c\t%d\t%d\n", last[i], W[i], Wm[i], colors[i]);

    fclose(boss_file);
}

double bwsd_expectation(int *t, int s, int n){
    int i;
	double value = 0.0;
    double frac;
    
	for(i = 1; i < n; i++){
        if(t[i] != 0){
            frac = (double)t[i]/s;
            value += i*frac;
        }
    }

    return value-1;
}

double bwsd_shannon_entropy(int *t, int s, int n){
    int i;
	double value = 0.0;
	
    for(i = 1; i < n; i++){
        if(t[i] != 0){
            double frac = (double)t[i]/(double)s;
            value += frac*log2(frac);
        }
    }

    return value*(-1);

}

void bwsd(int *DA, int n, double *expectation, double *entropy){
    int i;

    int *run_length = (int*)malloc((n*2)*sizeof(int));
    int current = 0;
    run_length[0] = 0;
    run_length[1] = 0;
    int pos = 1;
    for(i = 0; i < n; i++){
        if(DA[i] == current)
            run_length[pos]++;
        else {
            current = DA[i];
            run_length[pos+1]=current;
            run_length[pos+2]=1;
            pos += 2;
        }
    }
    if(run_length[pos-1] == 0){
        pos+= 2;
        run_length[pos-1] = 1;
        run_length[pos] = 0;
    }
    pos++; //size of run_lentgh

    int *t = (int*)malloc(n*sizeof(int));
    memset(t, 0, n*sizeof(int));
    for(i = 0; i < pos; i+=2){
        t[run_length[i+1]]++;
    }

    *expectation = bwsd_expectation(t, pos/2, n);
    *entropy = bwsd_shannon_entropy(t, pos/2, n);
}

void print_bwsd_matrixes(double **Dm, double **De, char **files, int files_n){
    int i,j;
    FILE *bwsd_matrixes = fopen("results/bwsd_matrixes.txt", "w");
    fprintf(bwsd_matrixes, "Expectation matrix (D_m):\n");

    fprintf(bwsd_matrixes, "\t");
    for(i = 0; i < files_n; i++)
        fprintf(bwsd_matrixes, "%10d\t", i+1);
    fprintf(bwsd_matrixes, "\n");

    for(i = 0; i < files_n; i++){
        fprintf(bwsd_matrixes, "%d\t", i+1);
        for(j = 0; j < files_n; j++){
            if(j >= i)
                fprintf(bwsd_matrixes, "%10.5lf\t", Dm[i][j]);
            else
                fprintf(bwsd_matrixes, "%10.5lf\t", Dm[j][i]);
        }
        fprintf(bwsd_matrixes, "\n");
    }
    fprintf(bwsd_matrixes, "\n\n");

    fprintf(bwsd_matrixes, "Shannon's entropy matrix (D_e):\n");

    fprintf(bwsd_matrixes, "\t");
    for(i = 0; i < files_n; i++)
        fprintf(bwsd_matrixes, "%10d\t", i+1);
    fprintf(bwsd_matrixes, "\n");

    for(i = 0; i < files_n; i++){
        fprintf(bwsd_matrixes, "%d\t", i+1);
        for(j = 0; j < files_n; j++){
            if(j >= i)
                fprintf(bwsd_matrixes, "%10.5lf\t", De[i][j]);
            else
                fprintf(bwsd_matrixes, "%10.5lf\t", De[j][i]);
        }
        fprintf(bwsd_matrixes, "\n");
    }
    
    fprintf(bwsd_matrixes, "\n\n");
    fprintf(bwsd_matrixes, "Matrix info:\n");
    for(i = 0; i < files_n; i++)
        fprintf(bwsd_matrixes, "%d: %s\n", i+1, files[i]);

    fprintf(bwsd_matrixes, "\n\n");
    fprintf(bwsd_matrixes, "CSV similarity matrix:\n");

    for(i = 0; i < files_n; i++)
        fprintf(bwsd_matrixes, ",%s",files[i]);
    fprintf(bwsd_matrixes, "\n");
    for(i = 0; i < files_n; i++){
        fprintf(bwsd_matrixes, "%s,", files[i]);
        for(j = 0; j < files_n; j++){
            if(j >= i)
                fprintf(bwsd_matrixes, "%lf,", Dm[i][j]);
            else
                fprintf(bwsd_matrixes, "%lf,", Dm[j][i]);
        }
        fprintf(bwsd_matrixes, "\n");
    }
    
    fprintf(bwsd_matrixes, "\n\n");
    fprintf(bwsd_matrixes, "CSV distance matrix:\n");

    for(i = 0; i < files_n; i++)
        fprintf(bwsd_matrixes, ",%s",files[i]);
    fprintf(bwsd_matrixes, "\n");
    for(i = 0; i < files_n; i++){
        fprintf(bwsd_matrixes, "%s,", files[i]);
        for(j = 0; j < files_n; j++){
            if(j >= i)
                fprintf(bwsd_matrixes, "%lf,", De[i][j]);
            else
                fprintf(bwsd_matrixes, "%lf,", De[j][i]);
        }
        fprintf(bwsd_matrixes, "\n");
    }
    

    fclose(bwsd_matrixes);
}