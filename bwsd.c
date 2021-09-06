#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <libgen.h>
#include "bwsd.h"

#define FILE_PATH 1024

#ifndef COVERAGE
	#define COVERAGE 0
#endif

//max length para o tamanho do arquivo

double bwsd_expectation(int *t, int s, int n){
    int i;
	double value = 0.0;
    double frac;
    
    //go to max_t
	for(i = 1; i < n; i++){
        if(t[i] != 0){
            frac = (double)t[i]/s;
            value += i*frac;
        }
    }

    return value-1.0;
}

double bwsd_shannon_entropy(int *t, int s, int n){
    int i;
	double value = 0.0;
	
    //go to max_t
    for(i = 1; i < n; i++){
        if(t[i] != 0){
            double frac = (double)t[i]/(double)s;
            value += frac*log2(frac);
        }
    }

    return value*(-1.0);

}

size_t apply_coverage_merge(int primaryCoverage, int secondaryCoverage, int *rl_freq, size_t pos){
    int repetitions = 0;
    while(repetitions <= secondaryCoverage){
        pos++;
        rl_freq[pos] = 1;
        repetitions++;
    }
    int remaining = primaryCoverage - repetitions;
    pos++;
    rl_freq[pos] = remaining;

    return pos;
}

void bwsd(char* file1, char* file2, size_t n, int k, double *expectation, double *entropy, int mem, int printBoss, char coverage_type){
    size_t i;

    size_t size = n+1;
    int *rl_freq = (int*)calloc(size, sizeof(int));
    int max_freq = 0;

    int current = 0;
    rl_freq[0] = 0;
    size_t pos = 0; // size of run_length

    short *colors = (short*)calloc((mem+2), sizeof(short));
    short *reduced_LCP = (short*)calloc((mem+2), sizeof(short));
    int *coverage = (int*)calloc((mem+2), sizeof(int));

    char color_file_name[128];
    char reduced_LCP_file_name[128];
    char coverage_file_name[128];

    sprintf(color_file_name, "results/%s-%s.2.colors", file1, file2);
    sprintf(reduced_LCP_file_name, "results/%s-%s.2.reduced_lcp", file1, file2);
    sprintf(coverage_file_name, "results/%s-%s.4.coverage", file1, file2);
    
    FILE *colors_file = fopen(color_file_name, "rb");
    FILE *recuced_LCP_file = fopen(reduced_LCP_file_name, "rb");
    FILE *coverage_file = fopen(coverage_file_name, "rb");

    fread(colors, sizeof(short), mem+1, colors_file);
    fread(reduced_LCP, sizeof(short), mem+1, colors_file);
    fread(coverage, sizeof(int), mem+1, coverage_file);

    int block_pos = 0;
    
    for(i = 0; i < n; i++){

        if(i != 0 && i%mem == 0){
            // todo: try to copy last value to position 0 and read mem elements from position 1 
            fseek(colors_file, -sizeof(short), SEEK_CUR);
            fread(colors, sizeof(short), mem+1, colors_file);

            fseek(recuced_LCP_file, -sizeof(short), SEEK_CUR);
            fread(reduced_LCP, sizeof(short), mem+1, recuced_LCP_file);

            fseek(coverage_file, -sizeof(int), SEEK_CUR);
            fread(coverage, sizeof(int), mem+1, coverage_file);

            block_pos=0;
        }

        #if COVERAGE 
        //if we have two same (k+1)-mers from distinct genomes, we break down their coverage frequencies and merge then intending to increase their similarity.

        //Ex 0^4 1^3 = 0^1 1^1 0^1 1^1 0^1 1^1 0^1
        if(coverage_type == 'e' && reduced_LCP[block_pos] >= k+1 && reduced_LCP[block_pos+1] >= k+1 && colors[block_pos] != colors[block_pos+1] && coverage[block_pos] > 1 ){
            if(coverage[block_pos] > coverage[block_pos+1]){
                pos = apply_coverage_merge(coverage[block_pos], coverage[block_pos+1], rl_freq, pos);
                if(rl_freq[pos] > max_freq)
                    max_freq = rl_freq[pos];
            } 
            else if(coverage[block_pos] < coverage[block_pos+1]) {
                pos = apply_coverage_merge(coverage[block_pos+1], coverage[block_pos], rl_freq, pos);
                if(rl_freq[pos] > max_freq)
                    max_freq = rl_freq[pos];
            } else {
                int repetitions = 0;
                while(repetitions < coverage[block_pos]){
                    pos++;
                    rl_freq[pos] = 1;
                    pos++;
                    rl_freq[pos] = 1;
                    repetitions++;
                }
            }
        } else if (coverage_type == 'd' && reduced_LCP[block_pos] != reduced_LCP[block_pos+1]) {
            if(colors[block_pos] == current)
                rl_freq[pos] += coverage[block_pos];
            else {
                if(rl_freq[pos] > max_freq)
                    max_freq = rl_freq[pos];

                current = colors[block_pos];
                pos++;
                rl_freq[pos]=1;
            }
        } else { 
        #endif
            if(colors[block_pos] == current){
                rl_freq[pos]++;
                #if COVERAGE
                    if(coverage_type == 'a' || (coverage_type == 'd' && reduced_LCP[block_pos] != reduced_LCP[block_pos+1]))
                        rl_freq[pos] += coverage[block_pos]-1;
                #endif
            } else {
                if(rl_freq[pos] > max_freq)
                    max_freq = rl_freq[pos];

                current = colors[block_pos];
                pos++;
                rl_freq[pos]=1;
            }
        #if COVERAGE
        }
        #endif

        block_pos++;
    }

    // if last color == 0, add 1^0
    if(colors[block_pos-1] == 0){
        pos++;
        rl_freq[pos] = 0;
    }

    // check if sum rl_freq = n;
    size_t sum_freq = 0;
    for(i = 0; i < pos; i++)
        sum_freq += rl_freq[i];
    // printf("diff n - sum_freq = %ld\n", n-sum_freq);

    free(colors); free(coverage); free(reduced_LCP);

    // computes every t_(k_j), where 1 <= j <= max_freq
    int *t = (int*) calloc((max_freq+1), sizeof(int));
    for(i = 1; i < pos; i++)
        t[rl_freq[i]]++;

    //sum termos

    free(rl_freq); 

    fclose(colors_file);
    fclose(recuced_LCP_file);
    fclose(coverage_file);

    char info[FILE_PATH];
    
    #if COVERAGE
        sprintf(info, "results/%s-%s_k_%d_coverage_1_%c.info", file1, file2, k, coverage_type);
    #else
        sprintf(info, "results/%s-%s_k_%d_coverage_0.info", file1, file2, k);
    #endif

    FILE *info_file = fopen(info, "w+");

    fprintf(info_file, "BWSD info of %s and %s genomes merge:\n", file1, file2);

    fprintf(info_file, "n = %ld\n", n);
    fprintf(info_file, "sum frequencies = %ld\n", sum_freq);
    fprintf(info_file, "n-sum_frequencies = %ld\n\n", n-sum_freq);

    fprintf(info_file, "s = %ld\n", pos);

    fprintf(info_file, "terms: \n");

    for(i = 0; i < max_freq+1; i++){
        if(t[i] != 0){
            fprintf(info_file, "t_%ld = %d\n", i, t[i]);
        }
    }

    fprintf(info_file, "\n");

    if(!printBoss){
        remove(color_file_name);
        remove(reduced_LCP_file_name);
        remove(coverage_file_name);
    }

    *expectation = bwsd_expectation(t, pos, max_freq);
    *entropy = bwsd_shannon_entropy(t, pos, max_freq);

    fprintf(info_file, "expectation = %lf\n", *expectation);
    fprintf(info_file, "entropy = %lf\n", *entropy);

    fclose(info_file);

    free(t);
}

void print_bwsd_matrixes(double **Dm, double **De, char **files, int files_n, char *path, int k, char coverage_type){
    int i,j;
    char *ptr;
    char outputFile[128];

    int len = strlen(path);
    char folder[len];
    strcpy(folder, basename(path));

    if(files_n > 2){
        ptr = strchr(path, '/');
        if (ptr != NULL)
            *ptr = '\0';

        #if COVERAGE
            sprintf(outputFile, "results/%s_distance_matrixes_k_%d_coverage_1_%c.txt", folder, k, coverage_type);
        #else
            sprintf(outputFile, "results/%s_distance_matrixes_k_%d_coverage_0.txt", folder, k);
        #endif
    } else {
        #if COVERAGE
            sprintf(outputFile, "results/%s-%s_distance_matrixes_k_%d_coverage_1_%c.txt", files[0], files[1], k, coverage_type);
        #else
            sprintf(outputFile, "results/%s-%s_distance_matrixes_k_%d_coverage_0.txt", files[0], files[1], k);
        #endif
    }
    

    FILE *bwsd_matrixes = fopen(outputFile, "w");
    
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
    fprintf(bwsd_matrixes, "CSV expectation distance matrix:\n");

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
    fprintf(bwsd_matrixes, "CSV entropy distance matrix:\n");

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

    fprintf(bwsd_matrixes, "\n\n");
    fprintf(bwsd_matrixes, "Phylip expectation distance matrix:\n");

    fprintf(bwsd_matrixes, "%d\n", files_n);
    for(i = 0; i < files_n; i++){
        fprintf(bwsd_matrixes, "%s", files[i]);
        for(j = 0; j < files_n; j++){
            if(j >= i)
                fprintf(bwsd_matrixes, "\t%lf", Dm[i][j]);
            else
                fprintf(bwsd_matrixes, "\t%lf", Dm[j][i]);
        }
        fprintf(bwsd_matrixes, "\n");
    }

    fprintf(bwsd_matrixes, "\n\n");
    fprintf(bwsd_matrixes, "Phylip entropy distance matrix:\n");

    fprintf(bwsd_matrixes, "%d\n", files_n);
    for(i = 0; i < files_n; i++){
        fprintf(bwsd_matrixes, "%s", files[i]);
        for(j = 0; j < files_n; j++){
            if(j >= i)
                fprintf(bwsd_matrixes, "\t%lf", De[i][j]);
            else
                fprintf(bwsd_matrixes, "\t%lf", De[j][i]);
        }
        fprintf(bwsd_matrixes, "\n");
    }

    fclose(bwsd_matrixes);
}
