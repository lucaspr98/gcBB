#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <libgen.h>
#include <time.h>
#include "bwsd.h"
#include "lib/rankbv.h"

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
    for(i = 1; i < n+1; i++){
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

    short *colors = (short*)calloc((mem+1), sizeof(short));
    short *summarized_LCP = (short*)calloc((mem+1), sizeof(short));
    short *summarized_SL = (short*)calloc((mem+1), sizeof(short));
    int *coverage = (int*)calloc((mem+1), sizeof(int));

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

    fread(colors, sizeof(short), mem, colors_file);
    fread(summarized_LCP, sizeof(short), mem, summarized_LCP_file);
    fread(summarized_SL, sizeof(short), mem, summarized_SL_file);
    fread(coverage, sizeof(int), mem, coverage_file);

    size_t *rl_freq = (size_t*)calloc(size, sizeof(size_t));
    size_t max_freq = 0;

    int current = consider1;
    rl_freq[consider1] = consider1;
    size_t pos = 0; // size of run_length
    int block_pos = 0;
    size_t rmq = 0;
    #if COVERAGE
    size_t consider1LastColorValue = consider1;
    size_t consider1LastCoverageValue = 0; 
    #endif

    for(i = 0; i < n; i++){     
        if(i != 0 && block_pos%mem == 0){
            fread(colors, sizeof(short), mem, colors_file);

            fread(summarized_LCP, sizeof(short), mem, summarized_LCP_file);

            fread(summarized_SL, sizeof(short), mem, summarized_SL_file);

            fread(coverage, sizeof(int), mem, coverage_file);

            block_pos=0;
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
        if(consider1LastColorValue == consider1 && colors[block_pos] == consider2 && rmq > k && (consider1LastCoverageValue > 1 || coverage[block_pos] > 1)){
            rl_freq[pos] = MAX((int)(rl_freq[pos])-1, 0); // decrease last 0 rl_freq because it is going to be intermixed with the current color
            pos++;
            rl_freq[pos++] = 0; // add 1^0 to rl_freq, since we are entering an intermix area and the last position is from genome 0
            apply_coverage_merge(consider1LastCoverageValue, coverage[block_pos], rl_freq, &pos);
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
                #if COVERAGE
                consider1LastColorValue = colors[block_pos];
                consider1LastCoverageValue = coverage[block_pos];
                #endif
            }
        #if COVERAGE
        }
        #endif
        block_pos++;
    }
    pos++;

    // check if sum rl_freq = n;
    // update max_freq;
    for(i = 0; i < pos; i++){
        max_freq = MAX(rl_freq[i], max_freq);
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
    fprintf(info_file, "pos = %ld\n\n", pos);

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

int getLastLCPGreaterThanKPos(short *lcp, int k, int intervalStart, int intervalEnd){
    int pos = intervalStart;
    for(int z = intervalStart+1; z < intervalEnd; z++){
        if(lcp[z] > k) {
            pos++;
        } else {
            break;
        }
    }
    return pos;
}


void bwsd_all(char* path, int samples, size_t n, size_t *sample_size, int k, int mem, double** Dm, double** De){
    size_t i, j, z;

    // Count computation time
    clock_t startClock, endClock;
    double cpu_time_used;

    startClock = clock();

    short *colors = (short*)calloc((mem+1), sizeof(short));
    short *summarized_LCP = (short*)calloc((mem+1), sizeof(short));
    short *summarized_SL = (short*)calloc((mem+1), sizeof(short));
    int *coverage = (int*)calloc((mem+1), sizeof(int));

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

    int tijSize = ((samples*(samples-1))/2)+1;

    size_t *lastJRank = calloc(tijSize, sizeof(size_t));
    size_t *lastIRank = calloc(tijSize, sizeof(size_t));

    size_t (**tij) = calloc(tijSize, sizeof(*tij));
    for(i = 0; i < samples-1; i++){
        for(j = i+1; j < samples; j++){
            int row = (((j-1)*(j))/2)+i;
            tij[row] = (size_t*) calloc(MAX(sample_size[i],sample_size[j])+2, sizeof(size_t));
        }
    } 
    size_t *tijMaxFreq = calloc(tijSize, sizeof(size_t));

    size_t *iCoverage = calloc(samples, sizeof(size_t));
    size_t *jCoverage = calloc(tijSize, sizeof(size_t));
    int needsToFindLcpNextBlock = 1;

    size_t blocks = ((n-1)/mem)+1;
    while(blocks){
        // last block
        int readSize = blocks == 1 && mem != n ? n%mem : mem; 
        fread(colors, sizeof(short), readSize, colors_file);
        fread(summarized_LCP, sizeof(short), readSize, summarized_LCP_file);
        fread(summarized_SL, sizeof(short), readSize, summarized_SL_file);
        fread(coverage, sizeof(int), readSize, coverage_file);

        rankbv_t **rbv = malloc(samples*sizeof(rankbv_t));
        for(i = 0; i < samples; i++){
            rbv[i] = rankbv_create(readSize, 2);
        }

        for(i = 0; i < readSize; i++){
            if(summarized_SL[i] > k) rankbv_setbit(rbv[colors[i]], i);
        }

        for(i = 0; i < samples; i++)
            rankbv_build(rbv[i]);

        for(i = 0; i < samples-1; i++){
            size_t intervalStart = 0;
            size_t intervalEnd = 0;

            for(z = 1; intervalEnd < readSize; z++){
                if(rankbv_access(rbv[i], intervalStart) == 1) iCoverage[i] = coverage[intervalStart];
                intervalEnd = rankbv_select1(rbv[i], z);
                // last interval of the block
                if(intervalEnd == -1) intervalEnd = readSize;
                 int lcpPos = -1;
                if(needsToFindLcpNextBlock){
                    lcpPos = getLastLCPGreaterThanKPos(summarized_LCP, k, intervalStart, intervalEnd);
                    if(lcpPos < intervalEnd && intervalEnd == readSize && rankbv_access(rbv[i], intervalEnd) == 1)
                        needsToFindLcpNextBlock = 0;
                }
                for(j = i+1; j < samples; j++){
                    int row = (((j-1)*(j))/2)+i;
                    size_t qtd;
                    int firstRbvJ1occurrence;
                    // if the following result is 0, we are in the
                    // start of a next block with unfinished interval
                    if(rankbv_access(rbv[i],intervalStart) == 1) 
                        firstRbvJ1occurrence = rankbv_select1(rbv[j], rankbv_rank1(rbv[j], intervalStart)+1);
                    else
                        firstRbvJ1occurrence = rankbv_select1(rbv[j], rankbv_rank1(rbv[j], intervalStart));
                    if(lcpPos >= intervalStart && firstRbvJ1occurrence >= intervalStart && firstRbvJ1occurrence <= lcpPos){
                        jCoverage[row] = coverage[firstRbvJ1occurrence];
                    }
                    // workaround for first interval, fail example:
                    // B_0 = 0 1 ...
                    // B_1 = 1 0 ...
                    if(intervalStart == 0) qtd = rankbv_rank1(rbv[j], intervalEnd);
                    else qtd = rankbv_rank1(rbv[j], intervalEnd)-rankbv_rank1(rbv[j], intervalStart);
                    // if we are looking the last interval of the block, 
                    // we store the qtd of the rbv[j]'s in lastJRank
                    if(intervalEnd == readSize && blocks != 1){
                        lastJRank[row] += qtd;
                        iCoverage[i] = coverage[intervalStart];
                    } else {
                        qtd += lastJRank[row];
                        lastJRank[row] = 0;
                        if(qtd != 0){
                            #if COVERAGE
                                if((jCoverage[row] > 0 && iCoverage[i] > 0) && (jCoverage[row] != 1 || iCoverage[i] != 1)){
                                    int commom =  MIN(jCoverage[row], iCoverage[i]);
                                    int difference = MAX(jCoverage[row], iCoverage[i]) - commom;
                                    tij[row][1] += commom*2;
                                    tij[row][difference]++;
                                    if(lastIRank[row] > 0) tij[row][lastIRank[row]-1]++;
                                    tij[row][qtd-1]++;
                                    tijMaxFreq[row] = MAX(tijMaxFreq[row], MAX(qtd-1,MAX(difference, lastIRank[row]-1)));
                                } else {
                            #endif 
                                tij[row][lastIRank[row]]++;
                                tij[row][qtd]++;
                                tijMaxFreq[row] = MAX(tijMaxFreq[row], MAX(qtd,lastIRank[row]));
                            #if COVERAGE
                                }
                            #endif
                            if(blocks == 1 && intervalEnd == readSize) 
                                lastIRank[row] = 0;
                            else 
                                lastIRank[row] = 1;
                        } else {
                            lastIRank[row]++;
                        }
                        jCoverage[row] = 0;
                        iCoverage[i] = 0;
                        needsToFindLcpNextBlock = 1;
                    }
                }
                intervalStart = intervalEnd;
            }
        }

        // update tij of lastIRank on last block
        if(blocks == 1){
            for(i = 0; i < samples-1; i++){
                for(j = i+1; j < samples; j++){
                    int row = (((j-1)*(j))/2)+i;
                    tij[row][lastIRank[row]]++;
                }
            }
        }

        for(i = 0; i < samples; i++) rankbv_free(rbv[i]);
        rankbv_free(rbv);

        blocks--;
    }

    free(lastJRank); free(lastIRank); free(iCoverage); free(jCoverage);

    char info[FILE_PATH];

    sprintf(info, "results/%s", path);

    #if COVERAGE
        strcat(info, "_coverage");
    #endif

    char extension[FILE_PATH];
    sprintf(extension, "_k_%d.info", k);
    strcat(info, extension);

    FILE *info_file = fopen(info, "a+");

    fprintf(info_file, "BWSD construction info of genomes from %s merge:\n\n", path);

    for(i = 0; i < samples-1; i++){
        for(j = i+1; j < samples; j++){
            int row = (((j-1)*(j))/2)+i;
            size_t s = 0;
            fprintf(info_file, "t_{%d,%d}\n", i,j);
            for(z = 1; z < tijMaxFreq[row]+1; z++){
                fprintf(info_file, "t_%d = %d\n", z, tij[row][z]);
                s += tij[row][z];
            }
            fprintf(info_file, "\n");
            Dm[j][i] = bwsd_expectation(tij[row], s, tijMaxFreq[row]);
            De[j][i] = bwsd_shannon_entropy(tij[row], s, tijMaxFreq[row]);
        }
    }

    for(i = 0; i < tijSize; i++)
        free(tij[i]);
    free(tij); 
    

    free(colors); free(coverage); free(summarized_LCP); free(summarized_SL); free(tijMaxFreq); 
    endClock = clock();

    cpu_time_used = ((double) (endClock - startClock)) / CLOCKS_PER_SEC;

    fprintf(info_file, "BWSD construction time: %lf seconds\n\n", cpu_time_used);

    fprintf(info_file, "\n");

    fclose(info_file);
}
