#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "boss.h"

#define FILE_PATH 1024

#ifndef COVERAGE
	#define COVERAGE 0
#endif

typedef struct {
    char W;
    short Wm, color;
    int coverage;
} kmer_range;

int compare(const void *element1, const void *element2) {
    kmer_range *e1 = (kmer_range *)element1; 
    kmer_range *e2 = (kmer_range *)element2;
    if(e1->W == e2->W) 
        return e1->color - e2->color;
    return e1->W - e2->W;
}

void Wi_sort(char *Wi, short *Wm, short *colors, int *coverage, int start, int end){
    int i;

    kmer_range *values = (kmer_range*)malloc(end*sizeof(kmer_range));
    for(i = start; i < end; i++){
        values[i].W = Wi[i];
        values[i].Wm = Wm[i];
        values[i].color = colors[i];
        values[i].coverage = coverage[i];
    }

    qsort(values, end, sizeof(kmer_range), compare);

    for(i = start; i < end; i++){
        Wi[i] = values[i].W;
        Wm[i] = values[i].Wm;
        colors[i] = values[i].color;
        coverage[i] = values[i].coverage;
    }
    
    free(values);
}

void add_edge(char *W, short **last, short *colors, short *summarized_LCP, short *summarized_SL, int freq, short *Wm, char bwt, int da, short lcp, short sl, int Wi_size, int edge_status){
    *W = bwt;
    *colors = da;
    *summarized_LCP = lcp;
    *summarized_SL = sl;
    if(edge_status == 0){
        if(Wi_size == 0){
            (*last)[Wi_size] = 1;
        } else {
            (*last)[Wi_size-1] = 0;
            (*last)[Wi_size] = 1;
        }
    } else if(edge_status == 1){
        (*last)[Wi_size] = 1;
    } else if(edge_status == 2){
        (*last)[Wi_size-1] = 0;
        (*last)[Wi_size] = 1;
    }

    if(freq == 0){
        *Wm = 1;
    }
}

size_t boss_construction(FILE *mergeLCP, FILE *mergeDA, FILE *mergeBWT, FILE *mergeSL, size_t n, int k, int samples, int mem, char* file1, char* file2, int printBoss, size_t *total_coverage){
    // Iterators
    size_t i = 0; // iterates through Wi
    int j = 0;
    size_t bi = 0; // iterates through BWT, LCP, SL and DA 
    int lcp_block_pos = 0; // iterates through LCP memory blocks
    int other_blocks_pos = 1; // iterates through BWT, SL and DA memory blocks

    // Count computation time
    clock_t start, end;
    double cpu_time_used;

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
    char boss_last[FILE_PATH];
    char boss_w[FILE_PATH];
    char boss_wm[FILE_PATH];
    char boss_colors[FILE_PATH];
    char boss_coverage[FILE_PATH];
    char boss_summarized_LCP[FILE_PATH];
    char boss_summarized_SL[FILE_PATH];

    sprintf(boss_last, "results/%s-%s_k_%d.2.last", file1, file2, k);
    sprintf(boss_w, "results/%s-%s_k_%d.1.W", file1, file2, k);
    sprintf(boss_wm, "results/%s-%s_k_%d.2.Wm", file1, file2, k);
    sprintf(boss_colors, "results/%s-%s_k_%d.2.colors", file1, file2, k);
    sprintf(boss_coverage, "results/%s-%s_k_%d.4.coverage", file1, file2, k);
    sprintf(boss_summarized_LCP, "results/%s-%s_k_%d.2.summarized_LCP", file1, file2, k);
    sprintf(boss_summarized_SL, "results/%s-%s_k_%d.2.summarized_SL", file1, file2, k);

    FILE *boss_last_file = fopen(boss_last, "wb");
    FILE *boss_w_file = fopen(boss_w, "wb");
    FILE *boss_wm_file = fopen(boss_wm, "wb");
    
    FILE *boss_colors_file = fopen(boss_colors, "wb");
    FILE *boss_coverage_file = fopen(boss_coverage, "wb");
    FILE *boss_summarized_LCP_file = fopen(boss_summarized_LCP, "wb");
    FILE *boss_summarized_SL_file = fopen(boss_summarized_SL, "wb");

    // BOSS construction variables
    short *last = (short*)calloc(50, sizeof(short));
    char *W = (char*)calloc(50, sizeof(char));
    short *Wm = (short*)calloc(50,sizeof(short));
    short *colors = (short*)calloc(50, sizeof(short));
    int *coverage = (int*)calloc(50, sizeof(int));
    short *summarized_LCP = (short*)calloc(50, sizeof(short));
    short *summarized_SL = (short*)calloc(50, sizeof(short));

    for(j = 0; j < 50; j++) coverage[j] = 1;

    int C[255] = { 0 };

    // BOSS construction auxiliary variables 
    int Wi_size = 0; 
    int W_freq[255] = { 0 }; // frequency of outgoing edges in a (k-1)-mer suffix range (detects W- = 1)
    int Wi_freq[255] = { 0 }; // frequency of outgoing edges in a k-mer suffix range (detects same outgoing edge in a vertex)
    int Wi_first_occurrence[samples][255]; // first occurence of an outgoing edge in a k-mer suffix range from a string collection
    memset(Wi_first_occurrence, 0,  sizeof(int)*samples*255);
    int DA_freq[samples][255]; // frequency of outgoing edges in a k-mer from a string collection (used to include same outgoing edge from distinct collections in BOSS representation)
    memset(DA_freq, 0,  sizeof(int)*samples*255);

    int dummies_freq[samples][255]; // frequency of outgoing edges from dummy inputs of size 1 ($)
    memset(dummies_freq, 0,  sizeof(int)*samples*255);

    while(bi < n){

        // read next block
        if(bi != 0 && lcp_block_pos%mem == 0){
            LCP[0] = LCP[mem];
            fread(LCP+1, sizeof(short), mem, mergeLCP);

            SL[0] = SL[mem-1]; SL[1] = SL[mem];
            fread(SL+2, sizeof(short), mem, mergeSL);
            
            DA[0] = DA[mem-1]; DA[1] = DA[mem];
            fread(DA+2, sizeof(char), mem, mergeDA);
            
            BWT[0] = BWT[mem-1]; BWT[1] = BWT[mem];
            fread(BWT+2, sizeof(char), mem, mergeBWT);
            
            for(j = 2; j < mem+3; j++) BWT[j] = (BWT[j] == 0) ? '$' : BWT[j];
            
            lcp_block_pos = 0;
            other_blocks_pos = 1;
        }

        // more than one outgoing edge of vertex i
        if(LCP[lcp_block_pos+1] >= k && bi != n-1 ){
            // since there is more than one outgoing edge, we don't need to check if BWT = $ or there is already BWT[bi] in Wi range
            // if(BWT[other_blocks_pos] != '$'){ //change BWT[other_blocks_pos] with a new last_BWT variable (?) removed for now
                if(Wi_freq[BWT[other_blocks_pos]] == 0){
                    // Add values to BOSS representation
                    add_edge(&W[Wi_size], &last, &colors[Wi_size], &summarized_LCP[Wi_size], &summarized_SL[Wi_size], W_freq[BWT[other_blocks_pos]], &Wm[Wi_size], BWT[other_blocks_pos], DA[other_blocks_pos], LCP[lcp_block_pos], SL[other_blocks_pos], Wi_size, 0);

                    Wi_first_occurrence[DA[other_blocks_pos]][BWT[other_blocks_pos]] = Wi_size;

                    // Increment variables
                    C[BWT[other_blocks_pos]]++; W_freq[BWT[other_blocks_pos]]++; Wi_freq[BWT[other_blocks_pos]]++; DA_freq[DA[other_blocks_pos]][BWT[other_blocks_pos]]++; Wi_size++; i++;
                    (*total_coverage)++;
                } else {
                    // check if there is already outgoing edge labeled with BWT[bi] from DA[bi] leaving vertex i
                    if(DA_freq[DA[other_blocks_pos]][BWT[other_blocks_pos]] == 0){
                        add_edge(&W[Wi_size], &last, &colors[Wi_size], &summarized_LCP[Wi_size], &summarized_SL[Wi_size], W_freq[BWT[other_blocks_pos]], &Wm[Wi_size], BWT[other_blocks_pos], DA[other_blocks_pos], LCP[lcp_block_pos], SL[other_blocks_pos], Wi_size, 0);

                        Wi_first_occurrence[DA[other_blocks_pos]][BWT[other_blocks_pos]] = Wi_size;

                        C[BWT[other_blocks_pos]]++; W_freq[BWT[other_blocks_pos]]++; Wi_freq[BWT[other_blocks_pos]]++; DA_freq[DA[other_blocks_pos]][BWT[other_blocks_pos]]++; Wi_size++; i++; 
                        (*total_coverage)++;
                    } else {
                        // increases the coverage information of the node with outgoing edge labeled with BWT[bi] from DA[bi] which is already on BOSS construction 
                        int existing_pos = Wi_first_occurrence[DA[other_blocks_pos]][BWT[other_blocks_pos]];
                        coverage[existing_pos]++;
                        (*total_coverage)++;
                    }
                }
            // }
        } else {
            // just one outgoing edge of vertex i
            if(Wi_size == 0){
                //fix SL[other_blocks_pos-1] memory leak
                if (SL[other_blocks_pos] == 1 && dummies_freq[DA[other_blocks_pos]][BWT[other_blocks_pos]] == 0) {
                    add_edge(&W[Wi_size], &last, &colors[Wi_size], &summarized_LCP[Wi_size], &summarized_SL[Wi_size], W_freq[BWT[other_blocks_pos]], &Wm[Wi_size], BWT[other_blocks_pos], DA[other_blocks_pos], LCP[lcp_block_pos], SL[other_blocks_pos], Wi_size, 1);

                    dummies_freq[DA[other_blocks_pos]][BWT[other_blocks_pos]]++;

                    C[BWT[other_blocks_pos]]++; W_freq[BWT[other_blocks_pos]]++; i++; Wi_size++;(*total_coverage)++;
                } else if(SL[other_blocks_pos] > 1 && !(LCP[lcp_block_pos] == SL[other_blocks_pos-1]-1 && BWT[other_blocks_pos] == BWT[other_blocks_pos-1] && DA[other_blocks_pos] == DA[other_blocks_pos-1])){
                    add_edge(&W[Wi_size], &last, &colors[Wi_size], &summarized_LCP[Wi_size], &summarized_SL[Wi_size], W_freq[BWT[other_blocks_pos]], &Wm[Wi_size], BWT[other_blocks_pos], DA[other_blocks_pos], LCP[lcp_block_pos], SL[other_blocks_pos], Wi_size, 1);
                    C[BWT[other_blocks_pos]]++; W_freq[BWT[other_blocks_pos]]++; i++; Wi_size++;(*total_coverage)++;
                } 
            } 
            // last outgoing edge of vertex i
            else {
                // check if there is already outgoing edge labeled with BWT[bi] leaving vertex i
                if(Wi_freq[BWT[other_blocks_pos]] == 0){
                    add_edge(&W[Wi_size], &last, &colors[Wi_size], &summarized_LCP[Wi_size], &summarized_SL[Wi_size], W_freq[BWT[other_blocks_pos]], &Wm[Wi_size], BWT[other_blocks_pos], DA[other_blocks_pos], LCP[lcp_block_pos], SL[other_blocks_pos], Wi_size, 2);

                    C[BWT[other_blocks_pos]]++; W_freq[BWT[other_blocks_pos]]++; Wi_size++; i++; 
                    (*total_coverage)++;
                } else {
                    // check if there is already outgoing edge labeled with BWT[bi] from DA[bi] leaving vertex i
                    if(DA_freq[DA[other_blocks_pos]][BWT[other_blocks_pos]] == 0){
                        add_edge(&W[Wi_size], &last, &colors[Wi_size], &summarized_LCP[Wi_size], &summarized_SL[Wi_size], W_freq[BWT[other_blocks_pos]], &Wm[Wi_size], BWT[other_blocks_pos], DA[other_blocks_pos], LCP[lcp_block_pos], SL[other_blocks_pos], Wi_size, 2);

                        C[BWT[other_blocks_pos]]++; W_freq[BWT[other_blocks_pos]]++; Wi_freq[BWT[other_blocks_pos]]++; DA_freq[DA[other_blocks_pos]][BWT[other_blocks_pos]]++; Wi_size++; i++;                   
                        (*total_coverage)++;
                    } else {
                        // increases the coverage information of the node with outgoing edge labeled with BWT[bi] from DA[bi] which is already on BOSS construction 
                        int existing_pos = Wi_first_occurrence[DA[other_blocks_pos]][BWT[other_blocks_pos]];
                        coverage[existing_pos]++;
                        (*total_coverage)++;
                    }
                }
                // sort outgoing edges of vertex i in lexigraphic order 
                if(Wi_size > 1){
                    Wi_sort(W, Wm, colors, coverage, 0, Wi_size);
                }                

                // clean frequency variables of outgoing edges in Wi 
                memset(Wi_freq, 0, sizeof(int)*255);   
                memset(DA_freq, 0, sizeof(int)*samples*255);
                memset(dummies_freq, 0, sizeof(int)*samples*255);
                memset(Wi_first_occurrence, 0, sizeof(int)*samples*255);
            }
            // if next LCP value is smaller than k-1 we have a new (k-1)-mer to keep track, so we clean W_freq values
            if(LCP[lcp_block_pos+1] < k-1){
                memset(W_freq, 0, sizeof(int)*255);
            }

            // Write Wi in BOSS results files
            if(printBoss){
                fwrite(last, sizeof(short), Wi_size, boss_last_file);
                fwrite(W, sizeof(char), Wi_size, boss_w_file);
                fwrite(Wm, sizeof(short), Wi_size, boss_wm_file);
            }

            // needed for bwsd computation
            fwrite(colors, sizeof(short), Wi_size, boss_colors_file);
            fwrite(coverage, sizeof(int), Wi_size, boss_coverage_file);
            fwrite(summarized_LCP, sizeof(short), Wi_size, boss_summarized_LCP_file);
            fwrite(summarized_SL, sizeof(short), Wi_size, boss_summarized_SL_file);

            // clean buffers
            memset(last, 0, sizeof(short)*50);   
            memset(W, 0, sizeof(char)*50);   
            memset(Wm, 0, sizeof(short)*50);   
            memset(colors, 0, sizeof(short)*50);   
            memset(coverage, 0, sizeof(int)*50);   
            memset(summarized_LCP, 0, sizeof(short)*50);   
            memset(summarized_SL, 0, sizeof(short)*50);   

            for(j = 0; j < 50; j++) coverage[j] = 1;

            Wi_size = 0; 
        }
        lcp_block_pos++;
        other_blocks_pos++;
        bi++;
    }

    // fix C values
    C[1] = C['$'];
    C[2] = C['A'] + C[1];
    C[3] = C['C'] + C[2];
    C[4] = C['G'] + C[3];
    C[5] = C['N'] + C[4];
    C[0] = 0;

    end = clock();

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    // Print BOSS construction info
    char alphabet[6] = {'$', 'A', 'C', 'G', 'N', 'T'};
    char info[FILE_PATH];

    sprintf(info, "results/%s-%s", file1, file2);

    #if COVERAGE
        strcat(info, "_coverage");
    #endif

    char extension[FILE_PATH];
    sprintf(extension, "_k_%d.info", k);
    strcat(info, extension);
    
    FILE *info_file = fopen(info, "w");
                
    fprintf(info_file, "Boss construction info of %s and %s genomes merge:\n\n", file1, file2);

    fprintf(info_file, "Boss construction time: %lf seconds\n\n", cpu_time_used);

    fprintf(info_file, "C array:\n");
    for(j = 0; j < 6; j++)
        fprintf(info_file, "%c %d\n", alphabet[j], C[j]);
    fprintf(info_file, "\n");

    fprintf(info_file, "Frequencies:\n");
    for(j = 0; j < 6; j++)
        fprintf(info_file, "%c %d\n", alphabet[j], C[alphabet[j]]);
    fprintf(info_file, "\n");

    fprintf(info_file, "Boss length: %ld\n\n", i);

    fprintf(info_file, "Total coverage: %ld\n\n", (*total_coverage));

    // Close used files
    fclose(info_file);
    
    if(printBoss){
        fclose(boss_last_file);
        fclose(boss_w_file);
        fclose(boss_wm_file);
    } else {
        remove(boss_last);
        remove(boss_w);
        remove(boss_wm);
    }

    fclose(boss_colors_file);
    fclose(boss_coverage_file);
    fclose(boss_summarized_LCP_file);
    fclose(boss_summarized_SL_file);

    // free BOSS construction needed variables
    free(LCP); free(BWT); free(DA); free(SL);
    
    // free BOSS construction variables
    free(last); free(W); free(Wm); free(colors); free(coverage); free(summarized_LCP); free(summarized_SL);

    return i;
};