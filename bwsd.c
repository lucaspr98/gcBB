#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <libgen.h>
#include <time.h>
#include "bwsd.h"

#define FILE_PATH 1024

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#ifndef COVERAGE
	#define COVERAGE 0
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

void apply_coverage_merge(int zeroCoverage, int oneCoverage, size_t *rl_freq, size_t *pos){
    while(zeroCoverage > 0 && oneCoverage > 0){
        rl_freq[(*pos)++] = 1;
        rl_freq[(*pos)++] = 1;
        zeroCoverage--;
        oneCoverage--;
    }
    int last = zeroCoverage == 0 ? 1 : 0;
    if(last == 1 && oneCoverage){
        rl_freq[(*pos)++] = 0;
        rl_freq[(*pos)++] = oneCoverage;
    } else if(zeroCoverage) {
        rl_freq[(*pos)++] = zeroCoverage;
        rl_freq[(*pos)++] = 0;
    }
    return;
}

void bwsd(char* path, size_t n, int k, double *expectation, double *entropy, int mem, size_t total_coverage, int consider1, int consider2){
    size_t i;

    #if COVERAGE
        size_t size = total_coverage+1;
    #else
       size_t size = n+1;
    #endif

    // Count computation time
    clock_t start, end;
    double cpu_time_used;

    start = clock();

    short *colors = (short*)calloc((mem+2), sizeof(short));
    short *summarized_LCP = (short*)calloc((mem+2), sizeof(short));
    short *summarized_SL = (short*)calloc((mem+2), sizeof(short));
    int *coverage = (int*)calloc((mem+2), sizeof(int));

    char color_file_name[FILE_PATH];
    char summarized_LCP_file_name[FILE_PATH];
    char summarized_SL_file_name[FILE_PATH];
    char coverage_file_name[FILE_PATH];

    sprintf(color_file_name, "results/%s_k_%d.2.colors", path, k);
    sprintf(summarized_LCP_file_name, "results/%s_k_%d.2.summarized_LCP", path, k);
    sprintf(summarized_SL_file_name, "results/%s_k_%d.2.summarized_SL", path, k);
    sprintf(coverage_file_name, "results/%s_k_%d.4.coverage", path, k);
    
    FILE *colors_file = fopen(color_file_name, "rb");
    FILE *summarized_LCP_file = fopen(summarized_LCP_file_name, "rb");
    FILE *summarized_SL_file = fopen(summarized_SL_file_name, "rb");
    FILE *coverage_file = fopen(coverage_file_name, "rb");

    fread(colors, sizeof(short), mem+1, colors_file);
    fread(summarized_LCP, sizeof(short), mem+1, summarized_LCP_file);
    fread(summarized_SL, sizeof(short), mem+1, summarized_SL_file);
    fread(coverage, sizeof(int), mem+1, coverage_file);

    size_t *rl_freq = (size_t*)calloc(size, sizeof(size_t));
    size_t max_freq = 0;

    int current = consider1;
    rl_freq[consider1] = consider1;
    size_t pos = 0; // size of run_length
    int block_pos = 0;
    
    for(i = 0; i < n; i++){     
        if(i != 0 && block_pos%mem == 0){
            colors[0] = colors[mem];
            fread(colors+1, sizeof(short), mem, colors_file);

            summarized_LCP[0] = summarized_LCP[mem];
            fread(summarized_LCP+1, sizeof(short), mem, summarized_LCP_file);

            summarized_SL[0] = summarized_SL[mem];
            fread(summarized_SL+1, sizeof(short), mem, summarized_SL_file);

            coverage[0] = coverage[mem];
            fread(coverage+1, sizeof(int), mem, coverage_file);

            block_pos=0;
        }

        if(colors[block_pos] != consider1 && colors[block_pos] != consider2) {
            block_pos++;
            continue;
        }

        if(colors[block_pos] != consider1 && colors[block_pos] != consider2) {
            rmq = MIN(rmq, summarized_LCP[block_pos]);
            block_pos++;
            continue;
        } else {
            rmq = summarized_LCP[block_pos];
        }

        #if COVERAGE 
        // If we have two same (k+1)-mers from distinct genomes, 
        // we break down their coverage frequencies and merge then 
        // intending to increase their similarity.

        // For instance, if a (k+1)-mer occurs 4 times in genome 0
        // and 3 times in genome 1 we have the following:
        // 0^4 1^3 = 0^1 1^1 0^1 1^1 0^1 1^1 0^1

        // We always intermix these value between 1^0 and 0^0 in 
        // order to "separate" the intermix from the "default" bwsd.
        // For example, 
        // ... 0^4 1^3 ... = ... 1^0 (0^1 1^1 0^1 1^1 0^1 1^1 0^1) 1^0 ...
        if(colors[block_pos-1] == 0 && colors[block_pos] == 1 && rmq > k && (coverage[block_pos-1] > 1 || coverage[block_pos] > 1)){
            rl_freq[pos++]--; // decrease last 0 rl_freq because it is going to be intermixed with the current color
            rl_freq[pos++] = 0; // add 1^0 to rl_freq, since we are entering an intermix area and the last position is from genome 0
            apply_coverage_merge(coverage[block_pos-1], coverage[block_pos], rl_freq, &pos);
            // set current to 0 to "restart" the bwsd 0s and 1s count
            current = 0;
        } else { 
        #endif
            if(summarized_SL[block_pos] > k){
                if(colors[block_pos] == current){
                    rl_freq[pos]++;
                } else {
                    current = colors[block_pos];
                    pos++;
                    rl_freq[pos]=1;
                }
            }
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
    // update max_freq;
    size_t sum_freq = 0;
    for(i = 0; i < pos; i++){
        max_freq = MAX(rl_freq[i], max_freq);
        sum_freq += rl_freq[i];
    }

    size_t s = pos;

    // update s removing 0's
    for(i = 0; i < pos; i++){
        if(rl_freq[i] == 0) s--;
    }

    free(colors); free(coverage); free(summarized_LCP); free(summarized_SL);

    // computes every t_(k_j), where 1 <= j <= max_freq
    size_t *t = (size_t*) calloc((max_freq+10), sizeof(size_t));
    short *genome0 = (short*) calloc((max_freq+10), sizeof(short));
    short *genome1 = (short*) calloc((max_freq+10), sizeof(short));
    for(i = 0; i < pos; i++){
        if(rl_freq[i] > 0){
            t[rl_freq[i]]++;
            if(i%2)
                genome0[rl_freq[i]] = 1;
            else 
                genome1[rl_freq[i]] = 1;
        } 
    }

    end = clock();

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    fclose(colors_file);
    fclose(summarized_LCP_file);
    fclose(summarized_SL_file);
    fclose(coverage_file);

    char info[FILE_PATH];

    sprintf(info, "results/%s", path);

    #if COVERAGE
        strcat(info, "_coverage");
    #endif

    char extension[FILE_PATH];
    sprintf(extension, "_k_%d.info", k);
    strcat(info, extension);

    FILE *info_file = fopen(info, "a+");

    fprintf(info_file, "BWSD info of %d and %d genomes merge:\n\n", consider1, consider2);

    fprintf(info_file, "BWSD construction time: %lf seconds\n\n", cpu_time_used);

    #if COVERAGE
        fprintf(info_file, "total_coverage = %ld\n", total_coverage);
    #else
        fprintf(info_file, "n = %ld\n", n);
    #endif
    fprintf(info_file, "sum frequencies = %ld\n", sum_freq);

    fprintf(info_file, "s = %ld\n\n", s);

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

    free(genome0); free(genome1);

    fprintf(info_file, "\n");
    free(rl_freq); 

    *expectation = bwsd_expectation(t, s, max_freq);
    *entropy = bwsd_shannon_entropy(t, s, max_freq);

    fprintf(info_file, "expectation = %lf\n", *expectation);
    fprintf(info_file, "entropy = %lf\n\n", *entropy);

    fclose(info_file);

    free(t);
}
