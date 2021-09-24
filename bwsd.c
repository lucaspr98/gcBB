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

#ifndef FILTER_CONTEXT
	#define FILTER_CONTEXT 0
#endif

double log2(double i){

	return log(i)/log(2);

}

double bwsd_expectation(size_t *t, size_t s, size_t n){
    size_t i;
	double value = 0.0;
    double frac;
    
    //n = max_t
	for(i = 1; i < n+1; i++){
        if(t[i] != 0){
            frac = (double)t[i]/s;
            value += i*frac;
        }
    }

    return value-1.0;
}

double bwsd_shannon_entropy(size_t *t, size_t s, size_t n){
    size_t i;
	double value = 0.0;
	
    //n = max_t
    for(i = 0; i < n+1; i++){
        if(t[i] != 0){
            double frac = (double)t[i]/(double)s;
            value += frac*(log2(frac));
        }
    }

    if(value)
        return value*(-1.0);
    return value;
}

size_t apply_coverage_merge(int primaryCoverage, int secondaryCoverage, size_t *rl_freq, size_t pos){
    int repetitions = 0;
    while(repetitions <= secondaryCoverage){
        pos++;
        rl_freq[pos] = 1;
        repetitions++;
    }
    int remaining = primaryCoverage - repetitions + 1;
    pos++;
    rl_freq[pos] = remaining;

    return pos;
}

void bwsd(char* file1, char* file2, size_t n, int k, double *expectation, double *entropy, int mem, int printBoss, char coverage_type, size_t total_coverage, short complement){
    size_t i, j;

    #if COVERAGE
        size_t size = total_coverage+1;
    #else
       size_t size = n+1;
    #endif

    short *colors = (short*)calloc((mem+2), sizeof(short));
    short *summarized_LCP = (short*)calloc((mem+2), sizeof(short));
    short *summarized_SL = (short*)calloc((mem+2), sizeof(short));
    int *coverage = (int*)calloc((mem+2), sizeof(int));

    char color_file_name[FILE_PATH];
    char summarized_LCP_file_name[FILE_PATH];
    char summarized_SL_file_name[FILE_PATH];
    char coverage_file_name[FILE_PATH];

    sprintf(color_file_name, "results/%s-%s.2.colors", file1, file2);
    sprintf(summarized_LCP_file_name, "results/%s-%s.2.summarized_LCP", file1, file2);
    sprintf(summarized_SL_file_name, "results/%s-%s.2.summarized_SL", file1, file2);
    sprintf(coverage_file_name, "results/%s-%s.4.coverage", file1, file2);
    
    FILE *colors_file = fopen(color_file_name, "rb");
    FILE *summarized_LCP_file = fopen(summarized_LCP_file_name, "rb");
    FILE *summarized_SL_file = fopen(summarized_SL_file_name, "rb");
    FILE *coverage_file = fopen(coverage_file_name, "rb");

    fread(colors, sizeof(short), mem+1, colors_file);
    if(complement) 
        for(j = 0; j < mem; j++) colors[j] = colors[j] == 1 ? 0 : 1;
    fread(summarized_LCP, sizeof(short), mem+1, summarized_LCP_file);
    fread(summarized_SL, sizeof(short), mem+1, summarized_SL_file);
    fread(coverage, sizeof(int), mem+1, coverage_file);

    size_t *rl_freq = (size_t*)calloc(size, sizeof(size_t));
    size_t max_freq = 0;

    int current = 0;
    rl_freq[0] = 0;
    size_t pos = 0; // size of run_length
    int block_pos = 0;
    
    for(i = 0; i < n; i++){

        if(i != 0 && i%mem == 0){
            // todo: try to copy last value to position 0 and read mem elements from position 1 
            fseek(colors_file, -sizeof(short), SEEK_CUR);
            fread(colors, sizeof(short), mem+1, colors_file);
            if(complement) 
                for(j = 0; j < mem; j++) colors[j] = colors[j] == 1 ? 0 : 1;

            fseek(summarized_LCP_file, -sizeof(short), SEEK_CUR);
            fread(summarized_LCP, sizeof(short), mem+1, summarized_LCP_file);

            fseek(summarized_SL_file, -sizeof(short), SEEK_CUR);
            fread(summarized_SL, sizeof(short), mem+1, summarized_SL_file);

            fseek(coverage_file, -sizeof(int), SEEK_CUR);
            fread(coverage, sizeof(int), mem+1, coverage_file);

            block_pos=0;
        }

        #if COVERAGE 
        //if we have two same (k+1)-mers from distinct genomes, we break down their coverage frequencies and merge then intending to increase their similarity.

        //Ex 0^4 1^3 = 0^1 1^1 0^1 1^1 0^1 1^1 0^1
        if(coverage_type == 'e' && summarized_LCP[block_pos] >= k+1 && summarized_LCP[block_pos+1] >= k+1 && colors[block_pos] != colors[block_pos+1] && coverage[block_pos] > 1 ){
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
        } else { 
        #endif
            #if FILTER_CONTEXT
            if(summarized_SL[block_pos] > k){
            #endif
                if(colors[block_pos] == current){
                    rl_freq[pos]++;
                    #if COVERAGE
                        if(coverage_type == 'a' || (coverage_type == 'd' && summarized_LCP[block_pos] != summarized_LCP[block_pos+1]))
                            rl_freq[pos] += coverage[block_pos]-1;
                    #endif
                } else {
                    if(rl_freq[pos] > max_freq)
                        max_freq = rl_freq[pos];

                    current = colors[block_pos];
                    pos++;
                    rl_freq[pos]=1;
                }
            #if FILTER_CONTEXT
            }
            #endif
        #if COVERAGE
        }
        #endif

        block_pos++;
    }
    pos++;

    // if last color == 0, add 1^0
    if(colors[block_pos-1] == 0){
        rl_freq[pos] = 0;
        pos++;
    } 

    // check if sum rl_freq = n;
    size_t sum_freq = 0;
    for(i = 0; i < pos; i++)
        sum_freq += rl_freq[i];

    free(colors); free(coverage); free(summarized_LCP); free(summarized_SL);

    // computes every t_(k_j), where 1 <= j <= max_freq
    size_t *t = (size_t*) calloc((max_freq+4), sizeof(size_t));
    short *genome0 = (short*) calloc((max_freq+4), sizeof(short));
    short *genome1 = (short*) calloc((max_freq+4), sizeof(short));
    for(i = 0; i < pos; i++){
        if(rl_freq[i] > 0){
            t[rl_freq[i]]++;
            if(i%2)
                genome0[rl_freq[i]] = 1;
            else 
                genome1[rl_freq[i]] = 1;
        } 
    }

    //sum termos

    fclose(colors_file);
    fclose(summarized_LCP_file);
    fclose(summarized_SL_file);
    fclose(coverage_file);

    char info[FILE_PATH];

    sprintf(info, "results/%s-%s_k_%d", file1, file2, k);

    #if COVERAGE
        char coverage_arg[FILE_PATH];
        sprintf(coverage_arg, "_coverage_1_%c", coverage_type);
        strcat(info, coverage_arg);
    #else
        strcat(info, "_coverage_0");
    #endif

    #if FILTER_CONTEXT
        strcat(info, "_filtered");
    #endif

    strcat(info, ".txt");

    FILE *info_file = fopen(info, "a+");

    fprintf(info_file, "BWSD info of %s and %s genomes merge:\n", file1, file2);

    #if COVERAGE
        fprintf(info_file, "total_coverage = %ld\n", total_coverage);
    #else
        fprintf(info_file, "n = %ld\n", n);
    #endif
    fprintf(info_file, "sum frequencies = %ld\n", sum_freq);
    #if COVERAGE
        fprintf(info_file, "total_coverage-sum_frequencies = %ld\n\n", total_coverage-sum_freq);
    #else   
        fprintf(info_file, "n-sum_frequencies = %ld\n\n", n-sum_freq);
    #endif

    fprintf(info_file, "s = %ld\n", pos);

    fprintf(info_file, "terms: \n");

    for(i = 0; i < max_freq+1; i++){
        if(t[i] != 0){
            fprintf(info_file, "t_%ld = %ld (", i, t[i]);
            if(genome0[i])
                fprintf(info_file, "0");
            if(genome0[i] && genome1[i])
                fprintf(info_file, ",");
            if(genome1[i])
                fprintf(info_file, "1");    
            fprintf(info_file, ")\n");
        }
    }

    fprintf(info_file, "\n");

    if(!printBoss){
        remove(color_file_name);
        remove(summarized_LCP_file_name);
        remove(summarized_SL_file_name);
        remove(coverage_file_name);
    }

    size_t s = pos;

    if(rl_freq[0] == 0) s--;
    if(rl_freq[pos-1] == 0) s--;

    free(rl_freq); 

    *expectation = bwsd_expectation(t, s, max_freq);
    *entropy = bwsd_shannon_entropy(t, s, max_freq);

    fprintf(info_file, "expectation = %lf\n", *expectation);
    fprintf(info_file, "entropy = %lf\n", *entropy);

    fclose(info_file);

    free(t);
}

void print_bwsd_matrixes(double **Dm, double **De, char **files, int files_n, char *path, int k, char coverage_type){
    int i,j;
    char *ptr;
    char outputFile[FILE_PATH];

    int len = strlen(path);
    char folder[len];
    strcpy(folder, basename(path));

    if(files_n > 2){
        ptr = strchr(path, '/');
        if (ptr != NULL)
            *ptr = '\0';

        sprintf(outputFile, "results/%s_distance_matrixes_k_%d", folder, k);

        #if COVERAGE
            char coverage_arg[FILE_PATH];
            sprintf(coverage_arg, "_coverage_1_%c", coverage_type);
            strcat(outputFile, coverage_arg);
        #else
            strcat(outputFile, "_coverage_0");
        #endif

        #if FILTER_CONTEXT
            strcat(outputFile, "_filtered");
        #endif

        strcat(outputFile, ".txt");

    } else {
        sprintf(outputFile, "results/%s-%s_distance_matrixes_k_%d", files[0], files[1], k);            

        char coverage_arg[FILE_PATH];
        #if COVERAGE
            sprintf(coverage_arg, "_coverage_1_%c", coverage_type);
            strcat(outputFile, coverage_arg);
        #else
            strcat(outputFile, "_coverage_0");
        #endif

        strcat(outputFile, ".txt");
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
            fprintf(bwsd_matrixes, "%10.5lf\t", Dm[i][j]);
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
            fprintf(bwsd_matrixes, "%10.5lf\t", De[i][j]);
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
            fprintf(bwsd_matrixes, "%lf,", Dm[i][j]);
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
            fprintf(bwsd_matrixes, "%lf,", De[i][j]);
        }
        fprintf(bwsd_matrixes, "\n");
    }

    fprintf(bwsd_matrixes, "\n\n");
    fprintf(bwsd_matrixes, "Phylip expectation distance matrix:\n");

    fprintf(bwsd_matrixes, "%d\n", files_n);
    for(i = 0; i < files_n; i++){
        fprintf(bwsd_matrixes, "%s", files[i]);
        for(j = 0; j < files_n; j++){
            fprintf(bwsd_matrixes, "\t%lf", Dm[i][j]);
        }
        fprintf(bwsd_matrixes, "\n");
    }

    fprintf(bwsd_matrixes, "\n\n");
    fprintf(bwsd_matrixes, "Phylip entropy distance matrix:\n");

    fprintf(bwsd_matrixes, "%d\n", files_n);
    for(i = 0; i < files_n; i++){
        fprintf(bwsd_matrixes, "%s", files[i]);
        for(j = 0; j < files_n; j++){
            fprintf(bwsd_matrixes, "\t%lf", De[i][j]);
        }
        fprintf(bwsd_matrixes, "\n");
    }

    fclose(bwsd_matrixes);
}
