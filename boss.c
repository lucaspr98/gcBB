#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "boss.h"
#include "external.h"

#define FILE_PATH 1024
#define ALPHABET_SIZE 255

typedef struct {
    char W;
    short Wm, color, summarizedLCP, summarizedSL;
    int coverage;
} kmerRange;

int compare(const void *element1, const void *element2) {
    kmerRange *e1 = (kmerRange *)element1; 
    kmerRange *e2 = (kmerRange *)element2;
    if(e1->W == e2->W) 
        return e1->color - e2->color;
    return e1->W - e2->W;
}

void WiSort(char *Wi, short *Wm, short *colors, int *coverage, short *summarizedSL, int start, int end){
    int i;

    kmerRange *values = (kmerRange*)malloc(end*sizeof(kmerRange));
    for(i = start; i < end; i++){
        values[i].W = Wi[i];
        values[i].Wm = Wm[i];
        values[i].color = colors[i];
        values[i].coverage = coverage[i];
        values[i].summarizedSL = summarizedSL[i];
    }

    qsort(values, end, sizeof(kmerRange), compare);

    for(i = start; i < end; i++){
        Wi[i] = values[i].W;
        Wm[i] = values[i].Wm;
        colors[i] = values[i].color;
        coverage[i] = values[i].coverage;
        summarizedSL[i] = values[i].summarizedSL;
    }
    
    free(values);
}

void fixWiLCP(char *W, short *summarizedLCP, int k, int WiSize){
    for(int i = 1; i < WiSize; i++){
        if(W[i-1] == W[i])
            summarizedLCP[i] = k+1;
        else 
            summarizedLCP[i] = k;
    }
}

void addEdge(char *W, short **last, short *colors, short *summarizedLCP, short *summarizedSL, int freq, short *Wm, char bwt, char da, short lcp, short sl, int WiSize, int edgeStatus){
    *W = bwt;
    *colors = da;
    *summarizedLCP = lcp;
    *summarizedSL = sl;
    if(edgeStatus == 0){
        if(WiSize == 0){
            (*last)[WiSize] = 1;
        } else {
            (*last)[WiSize-1] = 0;
            (*last)[WiSize] = 1;
        }
    } else if(edgeStatus == 1){
        (*last)[WiSize] = 1;
    } else if(edgeStatus == 2){
        (*last)[WiSize-1] = 0;
        (*last)[WiSize] = 1;
    }

    if(freq == 0){
        *Wm = 1;
    }
}

void bossConstruction(FILE *mergeLCP, FILE *mergeDA, FILE *mergeBWT, FILE *mergeSL, size_t n, int k, int samples, int mem, char* file1, char* file2, int printBoss){
    // Iterators
    unsigned long i = 0; // iterates through Wi
    int j = 0;
    size_t bi = 0; // iterates through BWT, LCP, SL and DA 
    int lcpBlockPos = 0; // iterates through LCP memory blocks
    int otherBlocksPos = 1; // iterates through BWT, SL and DA memory blocks

    // Count computation time
    clock_t start, end;
    double cpuTimeUsed;

    start = clock();

    // LCP, SL, DA and BWT blocks needed for BOSS construction
    short *LCP = (short*)calloc((mem+2), sizeof(short));
    short *SL = (short*)calloc((mem+3), sizeof(short));
    char *DA = (char*)calloc((mem+3), sizeof(char));
    char *BWT = (char*)calloc((mem+3), sizeof(char));

    fread(LCP, sizeof(short), mem+1, mergeLCP);
    fread(SL+1, sizeof(short), mem+1, mergeSL);
    fread(DA+1, sizeof(char), mem+1, mergeDA);
    fread(BWT+1, sizeof(char), mem+1, mergeBWT);
    for(j = 1; j < mem+3; j++) BWT[j] = (BWT[j] == 0) ? '$' : BWT[j];

    // BOSS result files
    char bossLast[FILE_PATH];
    char bossW[FILE_PATH];
    char bossWm[FILE_PATH];
    char bossColors[FILE_PATH];
    char bossCoverage[FILE_PATH];
    char bossSummarizedLCP[FILE_PATH];
    char bossSummarizedSL[FILE_PATH];

    #if ALL_VS_ALL
        snprintf(bossLast, FILE_PATH, "results/%s_k_%d.2.last", file1, k);
        snprintf(bossW, FILE_PATH, "results/%s_k_%d.1.W", file1, k);
        snprintf(bossWm, FILE_PATH, "results/%s_k_%d.2.Wm", file1, k);
        snprintf(bossColors, FILE_PATH, "results/%s_k_%d.2.colors", file1, k);
        snprintf(bossCoverage, FILE_PATH, "results/%s_k_%d.4.coverage", file1, k);
        snprintf(bossSummarizedLCP, FILE_PATH, "results/%s_k_%d.2.summarizedLCP", file1, k);
        snprintf(bossSummarizedSL, FILE_PATH, "results/%s_k_%d.2.summarizedSL", file1, k);
    #else
        snprintf(bossLast, FILE_PATH, "results/%s-%s_k_%d.2.last", file1, file2, k);
        snprintf(bossW, FILE_PATH, "results/%s-%s_k_%d.1.W", file1, file2, k);
        snprintf(bossWm, FILE_PATH, "results/%s-%s_k_%d.2.Wm", file1, file2, k);
        snprintf(bossColors, FILE_PATH, "results/%s-%s_k_%d.2.colors", file1, file2, k);
        snprintf(bossCoverage, FILE_PATH, "results/%s-%s_k_%d.4.coverage", file1, file2, k);
        snprintf(bossSummarizedLCP, FILE_PATH, "results/%s-%s_k_%d.2.summarizedLCP", file1, file2, k);
        snprintf(bossSummarizedSL, FILE_PATH, "results/%s-%s_k_%d.2.summarizedSL", file1, file2, k);
    #endif


    #if ALL_VS_ALL
        FILE *bossColorsFileExists = fopen(bossColors, "r");
        if(bossColorsFileExists){
            printf("BOSS needed files already computed\n");
            fclose(bossColorsFileExists);
            free(LCP); free(BWT); free(DA); free(SL);
            return;
        }
    #endif
    
    FILE *bossLastFile = fopen(bossLast, "wb");
    FILE *bossWFile = fopen(bossW, "wb");
    FILE *bossWm_file = fopen(bossWm, "wb");
    
    FILE *bossColorsFile = fopen(bossColors, "wb");
    FILE *bossCoverageFile = fopen(bossCoverage, "wb");
    FILE *bossSummarizedLCPFile = fopen(bossSummarizedLCP, "wb");
    FILE *bossSummarizedSLFile = fopen(bossSummarizedSL, "wb");

    // BOSS construction variables
    short *last = (short*)calloc(200, sizeof(short));
    char *W = (char*)calloc(200, sizeof(char));
    short *Wm = (short*)calloc(200,sizeof(short));
    short *colors = (short*)calloc(200, sizeof(short));
    int *coverage = (int*)calloc(200, sizeof(int));
    short *summarizedLCP = (short*)calloc(200, sizeof(short));
    short *summarizedSL = (short*)calloc(200, sizeof(short));

    for(j = 0; j < 200; j++) coverage[j] = 1;

    unsigned long C[ALPHABET_SIZE] = { 0 };

    // BOSS construction auxiliary variables 
    int WiSize = 0; 
    int WFreq[ALPHABET_SIZE] = { 0 }; // frequency of outgoing edges in a (k-1)-mer suffix range (detects W- = 1)
    int WiFreq[ALPHABET_SIZE] = { 0 }; // frequency of outgoing edges in a k-mer suffix range (detects same outgoing edge in a vertex)
    int WiFirstOccurrence[samples][ALPHABET_SIZE]; // first occurence of an outgoing edge in a k-mer suffix range from a string collection
    memset(WiFirstOccurrence, 0,  sizeof(int)*samples*ALPHABET_SIZE);
    int DAFreq[samples][ALPHABET_SIZE]; // frequency of outgoing edges in a k-mer from a string collection (used to include same outgoing edge from distinct collections in BOSS representation)
    memset(DAFreq, 0,  sizeof(int)*samples*ALPHABET_SIZE);

    int dummiesFreq[samples][ALPHABET_SIZE]; // frequency of outgoing edges from dummy inputs of size 1 ($)
    memset(dummiesFreq, 0,  sizeof(int)*samples*ALPHABET_SIZE);

    size_t *totalSampleColorsInBoss = calloc(samples, sizeof(size_t));
    size_t *totalSampleCoverageInBoss = calloc(samples, sizeof(size_t));

    while(bi < n){

        // read next block
        if(bi != 0 && lcpBlockPos%mem == 0){
            LCP[0] = LCP[mem];
            fread(LCP+1, sizeof(short), mem, mergeLCP);

            SL[0] = SL[mem]; SL[1] = SL[mem+1];
            fread(SL+2, sizeof(short), mem, mergeSL);
            
            DA[0] = DA[mem]; DA[1] = DA[mem+1];
            fread(DA+2, sizeof(char), mem, mergeDA);
            
            BWT[0] = BWT[mem]; BWT[1] = BWT[mem+1];
            fread(BWT+2, sizeof(char), mem, mergeBWT);
            
            for(j = 2; j < mem+3; j++) BWT[j] = (BWT[j] == 0) ? '$' : BWT[j];
            
            lcpBlockPos = 0;
            otherBlocksPos = 1;
        }

        // more than one outgoing edge of vertex i
        if(LCP[lcpBlockPos+1] >= k && bi != n-1 ){
            // since there is more than one outgoing edge, we don't need to check if BWT = $ or there is already BWT[bi] in Wi range
            if(WiFreq[BWT[otherBlocksPos]] == 0){
                // Add values to BOSS representation
                addEdge(&W[WiSize], &last, &colors[WiSize], &summarizedLCP[WiSize], &summarizedSL[WiSize], WFreq[BWT[otherBlocksPos]], &Wm[WiSize], BWT[otherBlocksPos], DA[otherBlocksPos], LCP[lcpBlockPos], SL[otherBlocksPos], WiSize, 0);
                WiFirstOccurrence[DA[otherBlocksPos]][BWT[otherBlocksPos]] = WiSize;
                // Increment variables
                C[BWT[otherBlocksPos]]++; WFreq[BWT[otherBlocksPos]]++; WiFreq[BWT[otherBlocksPos]]++; DAFreq[DA[otherBlocksPos]][BWT[otherBlocksPos]]++; WiSize++; i++;
                (totalSampleCoverageInBoss[DA[otherBlocksPos]])++;
                (totalSampleColorsInBoss[DA[otherBlocksPos]])++;
            } else {
                // check if there is already outgoing edge labeled with BWT[bi] from DA[bi] leaving vertex i
                if(DAFreq[DA[otherBlocksPos]][BWT[otherBlocksPos]] == 0){
                    addEdge(&W[WiSize], &last, &colors[WiSize], &summarizedLCP[WiSize], &summarizedSL[WiSize], WFreq[BWT[otherBlocksPos]], &Wm[WiSize], BWT[otherBlocksPos], DA[otherBlocksPos], LCP[lcpBlockPos], SL[otherBlocksPos], WiSize, 0);
                    WiFirstOccurrence[DA[otherBlocksPos]][BWT[otherBlocksPos]] = WiSize;
                    C[BWT[otherBlocksPos]]++; WFreq[BWT[otherBlocksPos]]++; WiFreq[BWT[otherBlocksPos]]++; DAFreq[DA[otherBlocksPos]][BWT[otherBlocksPos]]++; WiSize++; i++; 
                    (totalSampleCoverageInBoss[DA[otherBlocksPos]])++;
                    (totalSampleColorsInBoss[DA[otherBlocksPos]])++;
                } else {
                    // increases the coverage information of the node with outgoing edge labeled with BWT[bi] from DA[bi] which is already on BOSS construction 
                    int existingPos = WiFirstOccurrence[DA[otherBlocksPos]][BWT[otherBlocksPos]];
                    coverage[existingPos]++;
                    (totalSampleCoverageInBoss[DA[otherBlocksPos]])++;
                }
            }
        } else {
            // just one outgoing edge of vertex i
            if(WiSize == 0){
                //fix SL[otherBlocksPos-1] memory leak
                if (SL[otherBlocksPos] == 1 && dummiesFreq[DA[otherBlocksPos]][BWT[otherBlocksPos]] == 0) {
                    addEdge(&W[WiSize], &last, &colors[WiSize], &summarizedLCP[WiSize], &summarizedSL[WiSize], WFreq[BWT[otherBlocksPos]], &Wm[WiSize], BWT[otherBlocksPos], DA[otherBlocksPos], LCP[lcpBlockPos], SL[otherBlocksPos], WiSize, 1);

                    dummiesFreq[DA[otherBlocksPos]][BWT[otherBlocksPos]]++;

                    C[BWT[otherBlocksPos]]++; WFreq[BWT[otherBlocksPos]]++; i++; WiSize++;
                    
                    (totalSampleCoverageInBoss[DA[otherBlocksPos]])++;
                    (totalSampleColorsInBoss[DA[otherBlocksPos]])++;
                } else if(SL[otherBlocksPos] > 1 && !(LCP[lcpBlockPos] == SL[otherBlocksPos-1]-1 && BWT[otherBlocksPos] == BWT[otherBlocksPos-1] && DA[otherBlocksPos] == DA[otherBlocksPos-1])){
                    addEdge(&W[WiSize], &last, &colors[WiSize], &summarizedLCP[WiSize], &summarizedSL[WiSize], WFreq[BWT[otherBlocksPos]], &Wm[WiSize], BWT[otherBlocksPos], DA[otherBlocksPos], LCP[lcpBlockPos], SL[otherBlocksPos], WiSize, 1);
                    C[BWT[otherBlocksPos]]++; WFreq[BWT[otherBlocksPos]]++; i++; WiSize++;
                    
                    (totalSampleCoverageInBoss[DA[otherBlocksPos]])++;
                    (totalSampleColorsInBoss[DA[otherBlocksPos]])++;
                } 
            } 
            // last outgoing edge of vertex i
            else {
                // check if there is already outgoing edge labeled with BWT[bi] leaving vertex i
                if(WiFreq[BWT[otherBlocksPos]] == 0){
                    addEdge(&W[WiSize], &last, &colors[WiSize], &summarizedLCP[WiSize], &summarizedSL[WiSize], WFreq[BWT[otherBlocksPos]], &Wm[WiSize], BWT[otherBlocksPos], DA[otherBlocksPos], LCP[lcpBlockPos], SL[otherBlocksPos], WiSize, 2);

                    C[BWT[otherBlocksPos]]++; WFreq[BWT[otherBlocksPos]]++; WiSize++; i++; 
                    
                    (totalSampleCoverageInBoss[DA[otherBlocksPos]])++;
                    (totalSampleColorsInBoss[DA[otherBlocksPos]])++;
                } else {
                    // check if there is already outgoing edge labeled with BWT[bi] from DA[bi] leaving vertex i
                    if(DAFreq[DA[otherBlocksPos]][BWT[otherBlocksPos]] == 0){
                        addEdge(&W[WiSize], &last, &colors[WiSize], &summarizedLCP[WiSize], &summarizedSL[WiSize], WFreq[BWT[otherBlocksPos]], &Wm[WiSize], BWT[otherBlocksPos], DA[otherBlocksPos], LCP[lcpBlockPos], SL[otherBlocksPos], WiSize, 2);

                        C[BWT[otherBlocksPos]]++; WFreq[BWT[otherBlocksPos]]++; WiFreq[BWT[otherBlocksPos]]++; DAFreq[DA[otherBlocksPos]][BWT[otherBlocksPos]]++; WiSize++; i++;                   
                        
                        (totalSampleCoverageInBoss[DA[otherBlocksPos]])++;
                        (totalSampleColorsInBoss[DA[otherBlocksPos]])++;
                    } else {
                        // increases the coverage information of the node with outgoing edge labeled with BWT[bi] from DA[bi] which is already on BOSS construction 
                        int existingPos = WiFirstOccurrence[DA[otherBlocksPos]][BWT[otherBlocksPos]];
                        coverage[existingPos]++;
                        
                        (totalSampleCoverageInBoss[DA[otherBlocksPos]])++;
                    }
                }
                // sort outgoing edges of vertex i in lexigraphic order 
                if(WiSize > 1){
                    WiSort(W, Wm, colors, coverage, summarizedSL, 0, WiSize);
                    fixWiLCP(W, summarizedLCP, k, WiSize);
                }                

                // clean frequency variables of outgoing edges in Wi 
                memset(WiFreq, 0, sizeof(int)*ALPHABET_SIZE);   
                memset(DAFreq, 0, sizeof(int)*samples*ALPHABET_SIZE);
                memset(dummiesFreq, 0, sizeof(int)*samples*ALPHABET_SIZE);
                memset(WiFirstOccurrence, 0, sizeof(int)*samples*ALPHABET_SIZE);
            }
            // if next LCP value is smaller than k-1 we have a new (k-1)-mer to keep track, so we clean WFreq values
            if(LCP[lcpBlockPos+1] < k-1){
                memset(WFreq, 0, sizeof(int)*ALPHABET_SIZE);
            }

            // Write Wi in BOSS results files
            if(printBoss){
                fwrite(last, sizeof(short), WiSize, bossLastFile);
                fwrite(W, sizeof(char), WiSize, bossWFile);
                fwrite(Wm, sizeof(short), WiSize, bossWm_file);
            }

            // needed for bwsd computation
            fwrite(colors, sizeof(short), WiSize, bossColorsFile);
            fwrite(coverage, sizeof(int), WiSize, bossCoverageFile);
            fwrite(summarizedLCP, sizeof(short), WiSize, bossSummarizedLCPFile);
            fwrite(summarizedSL, sizeof(short), WiSize, bossSummarizedSLFile);

            // clean buffers
            memset(last, 0, sizeof(short)*200);   
            memset(W, 0, sizeof(char)*200);   
            memset(Wm, 0, sizeof(short)*200);   
            memset(colors, 0, sizeof(short)*200);   
            memset(coverage, 0, sizeof(int)*200);   
            memset(summarizedLCP, 0, sizeof(short)*200);   
            memset(summarizedSL, 0, sizeof(short)*200);   

            for(j = 0; j < 200; j++) coverage[j] = 1;

            WiSize = 0; 
        }
        lcpBlockPos++;
        otherBlocksPos++;
        bi++;
    }

    // fix C values
    C[1] = C['$'];
    C[2] = C['A'] + C[1];
    C[3] = C['C'] + C[2];
    C[4] = C['G'] + C[3];
    C[5] = C['N'] + C[4];
    C[0] = 0;

    #if ALL_VS_ALL
    FILE *bossInfoFile = getBossInfoFile(file1, NULL, k, 1);
    #else
    FILE *bossInfoFile = getBossInfoFile(file1, file2, k, 1);
    #endif

    fprintf(bossInfoFile, "%ld\n", i);
    for(j = 0; j < samples; j++){
        fprintf(bossInfoFile, "%ld ", totalSampleColorsInBoss[j]);
        // free(totalSampleColorsInBoss[j]);
    }
    free(totalSampleColorsInBoss);
    fprintf(bossInfoFile, "\n");
    for(j = 0; j < samples; j++){
        fprintf(bossInfoFile, "%ld ", totalSampleCoverageInBoss[j]);
        // free(totalSampleCoverageInBossm[j]);
    }
    free(totalSampleCoverageInBoss);
    fprintf(bossInfoFile, "\n");
    fclose(bossInfoFile);

    #if ALL_VS_ALL
    FILE *infoFile = getInfoFile(file1, NULL, k, 0);
    #else
    FILE *infoFile = getInfoFile(file1, file2, k, 0);
    #endif

    #if DEBUG
    char alphabet[6] = {'$', 'A', 'C', 'G', 'N', 'T'};
        #if ALL_VS_ALL
        printBOSSDebug(i, infoFile, file1, NULL, alphabet, C, totalSampleCoverageInBoss, samples);
        #else
        printBOSSDebug(i, infoFile, file1, file2, alphabet, C, totalSampleCoverageInBoss, samples);
        #endif
    #endif

    end = clock();

    cpuTimeUsed = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("BOSS construction time: %lf seconds\n", cpuTimeUsed);

    fprintf(infoFile, "BOSS construction time: %lf seconds\n", cpuTimeUsed);
    fclose(infoFile);

    // free BOSS construction needed variables
    free(LCP); free(BWT); free(DA); free(SL);
    
    // free BOSS construction variables
    free(last); free(W); free(Wm); free(colors); free(coverage); free(summarizedLCP); free(summarizedSL);
    
    if(printBoss){
        fclose(bossLastFile);
        fclose(bossWFile);
        fclose(bossWm_file);
    } else {
        fclose(bossLastFile);
        fclose(bossWFile);
        fclose(bossWm_file);
        remove(bossLast);
        remove(bossW);
        remove(bossWm);
    }

    fclose(bossColorsFile);
    fclose(bossCoverageFile);
    fclose(bossSummarizedLCPFile);
    fclose(bossSummarizedSLFile);

    return;
};

void printBOSSDebug(unsigned long bossLength, FILE* infoFile, char* file1, char* file2, char* alphabet, int* C, size_t* totalSampleCoverageInBoss, int samples){
    size_t j;
    #if ALL_VS_ALL
        fprintf(infoFile, "BOSS construction info of genomes from %s merge:\n\n", file1);
    #else
        fprintf(infoFile, "BOSS construction info of %s and %s genomes merge:\n\n", file1, file2);   
    #endif

    fprintf(infoFile, "C array:\n");
    for(j = 0; j < 6; j++)
        fprintf(infoFile, "%c %d\n", alphabet[j], C[j]);
    fprintf(infoFile, "\n");

    fprintf(infoFile, "Frequencies:\n");
    for(j = 0; j < 6; j++)
        fprintf(infoFile, "%c %d\n", alphabet[j], C[alphabet[j]]);
    fprintf(infoFile, "\n");

    fprintf(infoFile, "BOSS length: %ld\n\n", bossLength);

    unsigned long totalCoverage = 0;
    for(j = 0; j < samples; j++) totalCoverage += totalSampleCoverageInBoss[j];
    fprintf(infoFile, "Total coverage: %ld\n\n", totalCoverage);
}
