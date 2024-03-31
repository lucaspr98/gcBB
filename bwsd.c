#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <libgen.h>
#include <time.h>
#include "bwsd.h"
#include "external.h"
#include "lib/rankbv.h"

#define FILE_PATH 1024

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

typedef struct {
    unsigned long bossLen;
    size_t *totalSampleColorsInBoss;
    size_t *totalSampleCoverageInBoss;
} bossInfo;

double log2(double i){
	return log(i)/log(2);
}

double bwsdExpectation(size_t *t, size_t s, size_t n){
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

double bwsdShannonEntropy(size_t *t, size_t s, size_t n){
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

void applyCoverageMerge(int zeroCoverage, int oneCoverage, size_t *rlFreq, size_t *pos){
    while(zeroCoverage > 0 && oneCoverage > 0){
        rlFreq[(*pos)++] = 1;
        rlFreq[(*pos)++] = 1;
        zeroCoverage--;
        oneCoverage--;
    }
    int last = zeroCoverage == 0 ? 1 : 0;
    if(last == 1 && oneCoverage){
        rlFreq[(*pos)++] = 0;
        rlFreq[(*pos)++] = oneCoverage;
    } else if(zeroCoverage) {
        rlFreq[(*pos)++] = zeroCoverage;
        rlFreq[(*pos)++] = 0;
    }
    return;
}

bossInfo *getBossInfo(char* file1, char* file2, int k, int samples){
    FILE *infoFile = getBossInfoFile(file1, file2, k, 0);

    bossInfo *bossInfo = calloc(1,sizeof(*bossInfo));
    bossInfo->totalSampleColorsInBoss = calloc(samples, sizeof(size_t));
    bossInfo->totalSampleCoverageInBoss = calloc(samples, sizeof(size_t));

    fscanf(infoFile, "%ld", &(bossInfo->bossLen));
    for(int i = 0; i < samples; i++){
        fscanf(infoFile, "%ld", &(bossInfo->totalSampleColorsInBoss[i]));
    }
    for(int i = 0; i < samples; i++){
        fscanf(infoFile, "%ld", &(bossInfo->totalSampleCoverageInBoss[i]));
    }

    fclose(infoFile);

    return bossInfo;
}

void bwsd(char* file1, char* file2, int k, double *expectation, double *entropy, int mem, int printBoss, int consider1, int consider2){
    size_t i;

    bossInfo *info = getBossInfo(file1, file2, k, 2);
    unsigned long n = info->bossLen;
    size_t totalCoverage = info->totalSampleCoverageInBoss[consider1] +  info->totalSampleCoverageInBoss[consider2];
    #if COVERAGE
        size_t size = totalCoverage+1;
    #else
        size_t size = n+1;
    #endif

    // Count computation time
    clock_t start, end;
    double cpuTimeUsed;

    start = clock();

    short *colors = (short*)calloc((mem+1), sizeof(short));
    short *summarizedLCP = (short*)calloc((mem+1), sizeof(short));
    short *summarizedSL = (short*)calloc((mem+1), sizeof(short));
    int *coverage = (int*)calloc((mem+1), sizeof(int));

    char colorFileName[FILE_PATH];
    char summarizedLCPFileName[FILE_PATH];
    char summarizedSLFileName[FILE_PATH];
    char coverageFileName[FILE_PATH];

    sprintf(colorFileName, "results/%s-%s_k_%d.2.colors", file1, file2, k);
    sprintf(summarizedLCPFileName, "results/%s-%s_k_%d.2.summarizedLCP", file1, file2, k);
    sprintf(summarizedSLFileName, "results/%s-%s_k_%d.2.summarizedSL", file1, file2, k);
    sprintf(coverageFileName, "results/%s-%s_k_%d.4.coverage", file1, file2, k);
    
    FILE *colorsFile = fopen(colorFileName, "rb");
    FILE *summarizedLCPFile = fopen(summarizedLCPFileName, "rb");
    FILE *summarizedSLFile = fopen(summarizedSLFileName, "rb");
    FILE *coverageFile = fopen(coverageFileName, "rb");

    fread(colors, sizeof(short), mem, colorsFile);
    fread(summarizedLCP, sizeof(short), mem, summarizedLCPFile);
    fread(summarizedSL, sizeof(short), mem, summarizedSLFile);
    fread(coverage, sizeof(int), mem, coverageFile);

    size_t *rlFreq = (size_t*)calloc(size, sizeof(size_t));
    size_t maxFreq = 0;

    int current = consider1;
    rlFreq[consider1] = consider1;
    size_t pos = 0; // size of run_length
    int blockPos = 0;
    size_t rmq = 0;
    #if COVERAGE
    size_t consider1LastColorValue = consider1;
    size_t consider1LastCoverageValue = 0; 
    #endif


    for(int z= 0; z< n; z++){
        if(summarizedSL[z] > k)
            printf("%d ", colors[z]);
    }
    printf("\n");


    for(i = 0; i < n; i++){     
        if(i != 0 && blockPos%mem == 0){
            fread(colors, sizeof(short), mem, colorsFile);

            fread(summarizedLCP, sizeof(short), mem, summarizedLCPFile);

            fread(summarizedSL, sizeof(short), mem, summarizedSLFile);

            fread(coverage, sizeof(int), mem, coverageFile);

            blockPos=0;
        }

        if(colors[blockPos] != consider1 && colors[blockPos] != consider2) {
            rmq = MIN(rmq, summarizedLCP[blockPos]);
            blockPos++;
            continue;
        } else {
            rmq = summarizedLCP[blockPos];
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
        if(consider1LastColorValue == consider1 && colors[blockPos] == consider2 && rmq > k && (consider1LastCoverageValue > 1 || coverage[blockPos] > 1)){
            rlFreq[pos] = MAX((int)(rlFreq[pos])-1, 0); // decrease last 0 rlFreq because it is going to be intermixed with the current color
            pos++;
            rlFreq[pos++] = 0; // add 1^0 to rlFreq, since we are entering an intermix area and the last position is from genome 0
            applyCoverageMerge(consider1LastCoverageValue, coverage[blockPos], rlFreq, &pos);
            // set current to 0 to "restart" the bwsd 0s and 1s count
            current = 0;
        } else { 
        #endif
            if(summarizedSL[blockPos] > k){
                if(colors[blockPos] == current){
                    rlFreq[pos]++;
                } else {
                    current = colors[blockPos];
                    pos++;
                    rlFreq[pos]=1;
                }
                #if COVERAGE
                consider1LastColorValue = colors[blockPos];
                consider1LastCoverageValue = coverage[blockPos];
                #endif
            }
        #if COVERAGE
        }
        #endif
        blockPos++;
    }
    pos++;

    // check if sum rlFreq = n;
    // update maxFreq;
    for(i = 0; i < pos; i++){
        maxFreq = MAX(rlFreq[i], maxFreq);
    }

    size_t s = pos;

    // update s removing 0's
    for(i = 0; i < pos; i++){
        if(rlFreq[i] == 0) s--;
    }

    free(colors); free(coverage); free(summarizedLCP); free(summarizedSL);

    // computes every t_(k_j), where 1 <= j <= maxFreq
    size_t *t = (size_t*) calloc((maxFreq+10), sizeof(size_t));
    short *genome0 = (short*) calloc((maxFreq+10), sizeof(short));
    short *genome1 = (short*) calloc((maxFreq+10), sizeof(short));
    for(i = 0; i < pos; i++){
        if(rlFreq[i] > 0){
            t[rlFreq[i]]++;
            if(i%2)
                genome0[rlFreq[i]] = 1;
            else 
                genome1[rlFreq[i]] = 1;
        } 
    }

    *expectation = bwsdExpectation(t, s, maxFreq);
    *entropy = bwsdShannonEntropy(t, s, maxFreq);

    fclose(colorsFile);
    fclose(summarizedLCPFile);
    fclose(summarizedSLFile);
    fclose(coverageFile);
    
    if(!printBoss){
        // remove(colorFileName);
        // remove(summarizedLCPFileName);
        // remove(summarizedSLFileName);
        // remove(coverageFileName);
    }

    FILE* infoFile = getInfoFile(file1, file2, k, 1);

    #if DEBUG
    printBWSDDebug(infoFile, file1, file2, totalCoverage, n, pos, s, maxFreq, t, genome0, genome1);
    #endif

    free(rlFreq); 
    free(genome0); free(genome1);
    free(t);

    end = clock();

    cpuTimeUsed = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("BWSD computation time: %lf seconds\n", cpuTimeUsed);

    fprintf(infoFile, "BWSD computation time: %lf seconds\n", cpuTimeUsed);

    fclose(infoFile);

    return;
}

void printBWSDDebug(FILE* infoFile, char* file1, char* file2, size_t totalCoverage, size_t n, size_t pos, size_t s, size_t maxFreq, size_t* t, short* genome0, short* genome1){
    size_t i;

    fprintf(infoFile, "BWSD info of %s and %s genomes merge:\n\n", file1, file2);
    #if COVERAGE
        fprintf(infoFile, "totalCoverage = %ld\n", totalCoverage);
    #else
        fprintf(infoFile, "n = %ld\n", n);
    #endif
    fprintf(infoFile, "pos = %ld\n\n", pos);

    fprintf(infoFile, "s = %ld\n\n", s);

    fprintf(infoFile, "terms: \n");

    for(i = 0; i < maxFreq+1; i++){
        if(t[i] != 0){
            fprintf(infoFile, "t_%ld = %ld (", i, t[i]);
            if(genome0[i])
                fprintf(infoFile, "0");
            if(genome0[i] && genome1[i])
                fprintf(infoFile, ",");
            if(genome1[i])
                fprintf(infoFile, "1");    
            fprintf(infoFile, ")\n");
        }
    }
    fprintf(infoFile, "\n");
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

void bwsdAll(char* path, int samples, int k, int mem, double** Dm, double** De){
    size_t i, j, z;

    // Count computation time
    clock_t startClock, endClock;
    double cpuTimeUsed;

    startClock = clock();

    bossInfo *info = getBossInfo(path, NULL, k, samples);
    unsigned long n = info->bossLen;
    #if COVERAGE
        size_t *sampleSize = info->totalSampleCoverageInBoss;
    #else
        size_t *sampleSize = info->totalSampleColorsInBoss;
    #endif
    
    short *colors = (short*)calloc((mem+1), sizeof(short));
    short *summarizedLCP = (short*)calloc((mem+1), sizeof(short));
    short *summarizedSL = (short*)calloc((mem+1), sizeof(short));
    int *coverage = (int*)calloc((mem+1), sizeof(int));

    char colorFileName[FILE_PATH];
    char summarizedLCPFileName[FILE_PATH];
    char summarizedSLFileName[FILE_PATH];
    char coverageFileName[FILE_PATH];

    sprintf(colorFileName, "results/%s_k_%d.2.colors", path, k);
    sprintf(summarizedLCPFileName, "results/%s_k_%d.2.summarizedLCP", path, k);
    sprintf(summarizedSLFileName, "results/%s_k_%d.2.summarizedSL", path, k);
    sprintf(coverageFileName, "results/%s_k_%d.4.coverage", path, k);
    
    FILE *colorsFile = fopen(colorFileName, "rb");
    FILE *summarizedLCPFile = fopen(summarizedLCPFileName, "rb");
    FILE *summarizedSLFile = fopen(summarizedSLFileName, "rb");
    FILE *coverageFile = fopen(coverageFileName, "rb");

    int tijSize = ((samples*(samples-1))/2)+1;

    size_t *lastJRank = calloc(tijSize, sizeof(size_t));
    size_t *lastIRank = calloc(tijSize, sizeof(size_t));

    size_t (**tij) = calloc(tijSize, sizeof(*tij));
    for(i = 0; i < samples-1; i++){
        for(j = i+1; j < samples; j++){
            int row = (((j-1)*(j))/2)+i;
            tij[row] = (size_t*) calloc(MAX(sampleSize[i],sampleSize[j])+2, sizeof(size_t));
        }
    } 
    size_t *tijMaxFreq = calloc(tijSize, sizeof(size_t));

    size_t *iCoverage = calloc(samples, sizeof(size_t));
    size_t *jCoverage = calloc(tijSize+1, sizeof(size_t));
    int needsToFindLcpNextBlock = 1;

    size_t blocks = ((n-1)/mem)+1;

    while(blocks){
        // last block
        int readSize = blocks == 1 && mem != n ? n%mem : mem; 
        fread(colors, sizeof(short), readSize, colorsFile);
        fread(summarizedLCP, sizeof(short), readSize, summarizedLCPFile);
        fread(summarizedSL, sizeof(short), readSize, summarizedSLFile);
        fread(coverage, sizeof(int), readSize, coverageFile);

        for(int z = 0; z < samples-1; z++){
            for(int x = z+1; x < samples; x++){
                for(int y = 0; y < n; y++){
                    if((colors[y] == z || colors[y] == x) && summarizedSL[y] > k)
                        printf("%d ", colors[y]);
                }
                printf("\n");
                for(int y = 0; y < n; y++){
                    if((colors[y] == z || colors[y] == x) && summarizedSL[y] > k)
                        printf("%d ", summarizedSL[y]);
                }
                printf("\n");

                for(int y = 0; y < n; y++){
                    if((colors[y] == z || colors[y] == x) && summarizedSL[y] > k)
                        printf("%d ", summarizedLCP[y]);
                }
                printf("\n");
            }
        }

        rankbv_t **rbv = malloc(samples*sizeof(rankbv_t));
        for(i = 0; i < samples; i++){
            rbv[i] = rankbv_create(readSize, 2);
        }

        for(i = 0; i < readSize; i++){
            if(summarizedSL[i] > k) {
                rankbv_setbit(rbv[colors[i]], i);
            }
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
                    lcpPos = getLastLCPGreaterThanKPos(summarizedLCP, k, intervalStart, intervalEnd);
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
                    if(lcpPos >= intervalStart && firstRbvJ1occurrence >= intervalStart && firstRbvJ1occurrence <= lcpPos && firstRbvJ1occurrence != -1){
                        jCoverage[row] = coverage[firstRbvJ1occurrence];
                    } else {
                        jCoverage[row] = 0;
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
        free(rbv);

        blocks--;
    }

    for(i = 0; i < samples-1; i++){
        for(j = i+1; j < samples; j++){
            int row = (((j-1)*(j))/2)+i;
            size_t s = 0;
            for(z = 1; z < tijMaxFreq[row]+1; z++) s += tij[row][z];
            Dm[j][i] = bwsdExpectation(tij[row], s, tijMaxFreq[row]);
            De[j][i] = bwsdShannonEntropy(tij[row], s, tijMaxFreq[row]);
        }
    }

    FILE *infoFile = getInfoFile(path, NULL, k, 1);

    #if DEBUG
    printBWSDALLDebug(infoFile, path, samples, tijMaxFreq, tij);
    #endif

    endClock = clock();

    cpuTimeUsed = ((double) (endClock - startClock)) / CLOCKS_PER_SEC;

    printf("BWSD computation time: %lf seconds\n", cpuTimeUsed);

    fprintf(infoFile, "BWSD computation time: %lf seconds\n", cpuTimeUsed);
    fclose(infoFile);

    free(lastJRank); free(lastIRank); free(iCoverage); free(jCoverage);

    for(i = 0; i < tijSize; i++)
        free(tij[i]);
    free(tij); 

    free(colors); free(coverage); free(summarizedLCP); free(summarizedSL); free(tijMaxFreq); 

    free(info->totalSampleColorsInBoss);
    free(info->totalSampleCoverageInBoss);
    free(info);

    fclose(colorsFile);
    fclose(summarizedLCPFile);
    fclose(summarizedSLFile);
    fclose(coverageFile);

    

    return;
}

void printBWSDALLDebug(FILE* infoFile, char* path, int samples, size_t* tijMaxFreq,size_t** tij){
    size_t i, j, z;
    fprintf(infoFile, "BWSD computation info of genomes from %s merge:\n\n", path);    
    for(i = 0; i < samples-1; i++){
        for(j = i+1; j < samples; j++){
            int row = (((j-1)*(j))/2)+i;
            fprintf(infoFile, "t_{%ld,%ld}\n", i,j);
            for(z = 1; z < tijMaxFreq[row]+1; z++){
                if(tij[row][z] != 0){
                    fprintf(infoFile, "t_%ld = %ld\n", z, tij[row][z]);
                }
            }
            fprintf(infoFile, "\n");
        }
    }
    fprintf(infoFile, "\n");
}
