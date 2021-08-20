#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "boss.h"

#define FILE_PATH 1024

void Wi_sort(char *Wi, short *Wm, short *colors, int *coverage, int start, int end){
    int i;
    int range = end-start;
    char Wi_tmp[range];
    short Wm_tmp[range];
    short colors_tmp[range];
    int coverage_tmp[range];
    int Wi_aux[255];
    int Wm_aux[255];
    int repetitive[255];

    memset(Wi_tmp, 0, sizeof(char)*(range));    
    memset(Wm_tmp, 0, sizeof(int)*(range)); 
    memset(colors_tmp, 0, sizeof(int)*(range)); 
    memset(coverage_tmp, 0, sizeof(int)*(range)); 
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
        coverage_tmp[Wi_aux[Wi[i]]-1] = coverage[i];
        Wi_aux[Wi[i]]--;
    }

    for(i = start; i < end; i++){
        Wi[i] = Wi_tmp[i-start];
        Wm[i] = Wm_tmp[i-start];
        colors[i] = colors_tmp[i-start];
        coverage[i] = coverage_tmp[i-start];
    }
}

void add_edge(char *W, short **last, short *colors, short *reduced_LCP, int freq, short *Wm, char bwt, int da, short lcp, int Wi_size, int edge_status){
    *W = bwt;
    *colors = da;
    *reduced_LCP = lcp;
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

size_t boss_construction(FILE *mergeLCP, FILE *mergeDA, FILE *mergeBWT, size_t n, int k, int samples, int mem, char* file1, char* file2, int printBoss){
    // Iterators
    size_t i = 0; // iterates through Wi
    int j = 0;
    size_t bi = 0; // iterates through BWT, LCP and DA 
    int block_pos = 0; // iterates through BWT, LCP and DA memory blocks

    // LCP, DA and BWT blocks needed for BOSS construction
    short *LCP = (short*)calloc((mem+2), sizeof(short));
    char *DA = (char*)calloc((mem+1), sizeof(char));
    char *BWT = (char*)calloc((mem+1), sizeof(char));

    fread(LCP, sizeof(short), mem+1, mergeLCP);
    fread(DA, sizeof(char), mem, mergeDA);
    fread(BWT, sizeof(char), mem, mergeBWT);
    for(j = 0; j < mem; j++) BWT[j] = (BWT[j] == 0) ? '$' : BWT[j];

    // BOSS result files
    char boss_last[FILE_PATH];
    char boss_w[FILE_PATH];
    char boss_wm[FILE_PATH];
    char boss_colors[FILE_PATH];
    char boss_coverage[FILE_PATH];
    char boss_reduced_LCP[FILE_PATH];

    sprintf(boss_last, "results/%s-%s.2.last", file1, file2);
    sprintf(boss_w, "results/%s-%s.1.W", file1, file2);
    sprintf(boss_wm, "results/%s-%s.2.Wm", file1, file2);
    sprintf(boss_colors, "results/%s-%s.2.colors", file1, file2);
    sprintf(boss_coverage, "results/%s-%s.4.coverage", file1, file2);
    sprintf(boss_reduced_LCP, "results/%s-%s.2.reduced_lcp", file1, file2);

    FILE *boss_last_file = fopen(boss_last, "wb");
    FILE *boss_w_file = fopen(boss_w, "wb");
    FILE *boss_wm_file = fopen(boss_wm, "wb");
    
    FILE *boss_colors_file = fopen(boss_colors, "wb");
    FILE *boss_coverage_file = fopen(boss_coverage, "wb");
    FILE *boss_reduced_LCP_file = fopen(boss_reduced_LCP, "wb");

    // BOSS construction variables
    short *last = (short*)calloc(10, sizeof(short));
    char *W = (char*)calloc(10, sizeof(char));
    short *Wm = (short*)calloc(10,sizeof(short));
    short *colors = (short*)calloc(10, sizeof(short));
    int *coverage = (int*)calloc(10, sizeof(int));
    short *reduced_LCP = (short*)calloc(10, sizeof(short));

    for(j = 0; j < 10; j++) coverage[j] = 1;

    int C[255] = { 0 };
    size_t total_coverage = 0;

    // BOSS construction auxiliary variables 
    int Wi_size = 0; 
    int W_freq[255] = { 0 }; // frequency of outgoing edges in a (k-1)-mer suffix range (detects W- = 1)
    int Wi_freq[255] = { 0 }; // frequency of outgoing edges in a k-mer suffix range (detects same outgoing edge in a vertex)
    int Wi_first_occurrence[samples][255]; // first occurence of an outgoing edge in a k-mer suffix range from a string collection
    memset(Wi_first_occurrence, 0,  sizeof(int)*samples*255);
    int DA_freq[samples][255]; // frequency of outgoing edges in a k-mer from a string collection (used to include same outgoing edge from distinct collections in BOSS representation)
    memset(DA_freq, 0,  sizeof(int)*samples*255);

    while(bi < n){

        // read next block
        if(bi != 0 && bi%mem == 0){
            fseek(mergeLCP, -2L, SEEK_CUR);
            fread(LCP, sizeof(short), mem+1, mergeLCP);
            fread(DA, sizeof(char), mem, mergeDA);
            fread(BWT, sizeof(char), mem, mergeBWT);
            for(j = 0; j < mem; j++) BWT[j] = (BWT[j] == 0) ? '$' : BWT[j];
            
            block_pos = 0;
        }

        // more than one outgoing edge of vertex i
        if(LCP[block_pos+1] >= k && bi != n-1){
            // since there is more than one outgoing edge, we don't need to check if BWT = $ or there is already BWT[bi] in Wi range
            if(BWT[block_pos] != '$'){ //change BWT[block_pos] with a new last_BWT variable
                if(Wi_freq[BWT[block_pos]] == 0){
                    // Add values to BOSS representation
                    add_edge(&W[Wi_size], &last, &colors[Wi_size], &reduced_LCP[Wi_size], W_freq[BWT[block_pos]], &Wm[Wi_size], BWT[block_pos], DA[block_pos], LCP[block_pos], Wi_size, 0);

                    Wi_first_occurrence[DA[block_pos]][BWT[block_pos]] = Wi_size;

                    // Increment variables
                    C[BWT[block_pos]]++; W_freq[BWT[block_pos]]++; Wi_freq[BWT[block_pos]]++; DA_freq[DA[block_pos]][BWT[block_pos]]++; Wi_size++; i++;
                    total_coverage++;
                } else {
                    // check if there is already outgoing edge labeled with BWT[bi] from DA[bi] leaving vertex i
                    if(DA_freq[DA[block_pos]][BWT[block_pos]] == 0){
                        add_edge(&W[Wi_size], &last, &colors[Wi_size], &reduced_LCP[Wi_size], W_freq[BWT[block_pos]], &Wm[Wi_size], BWT[block_pos], DA[block_pos], LCP[block_pos], Wi_size, 0);

                        Wi_first_occurrence[DA[block_pos]][BWT[block_pos]] = Wi_size;

                        C[BWT[block_pos]]++; W_freq[BWT[block_pos]]++; Wi_freq[BWT[block_pos]]++; DA_freq[DA[block_pos]][BWT[block_pos]]++; Wi_size++; i++; 
                        total_coverage++;
                    } else {
                        // increases the coverage information of the node with outgoing edge labeled with BWT[bi] from DA[bi] which is already on BOSS construction 
                        int existing_pos = Wi_first_occurrence[DA[block_pos]][BWT[block_pos]];
                        coverage[existing_pos]++;
                        total_coverage++;
                    }
                }
            }
        } else {
            // just one outgoing edge of vertex i
            if(Wi_size == 0){
                add_edge(&W[Wi_size], &last, &colors[Wi_size], &reduced_LCP[Wi_size], W_freq[BWT[block_pos]], &Wm[Wi_size], BWT[block_pos], DA[block_pos], LCP[block_pos], Wi_size, 1);

                C[BWT[block_pos]]++; W_freq[BWT[block_pos]]++; i++; Wi_size++;
                total_coverage++;
            } 
            // last outgoing edge of vertex i
            else {
                // check if there is already outgoing edge labeled with BWT[bi] leaving vertex i
                if(Wi_freq[BWT[block_pos]] == 0){
                    add_edge(&W[Wi_size], &last, &colors[Wi_size], &reduced_LCP[Wi_size], W_freq[BWT[block_pos]], &Wm[Wi_size], BWT[block_pos], DA[block_pos], LCP[block_pos], Wi_size, 2);

                    C[BWT[block_pos]]++; W_freq[BWT[block_pos]]++; Wi_size++; i++; 
                    total_coverage++;
                } else {
                    // check if there is already outgoing edge labeled with BWT[bi] from DA[bi] leaving vertex i
                    if(DA_freq[DA[block_pos]][BWT[block_pos]] == 0){
                        add_edge(&W[Wi_size], &last, &colors[Wi_size], &reduced_LCP[Wi_size], W_freq[BWT[block_pos]], &Wm[Wi_size], BWT[block_pos], DA[block_pos], LCP[block_pos], Wi_size, 2);

                        C[BWT[block_pos]]++; W_freq[BWT[block_pos]]++; Wi_freq[BWT[block_pos]]++; DA_freq[DA[block_pos]][BWT[block_pos]]++; Wi_size++; i++;                   
                        total_coverage++;
                    } else {
                        // increases the coverage information of the node with outgoing edge labeled with BWT[bi] from DA[bi] which is already on BOSS construction 
                        int existing_pos = Wi_first_occurrence[DA[block_pos]][BWT[block_pos]];
                        coverage[existing_pos]++;
                        total_coverage++;
                    }
                }
                // sort outgoing edges of vertex i in lexigraphic order 
                if(Wi_size > 1){
                    Wi_sort(W, Wm, colors, coverage, 0, Wi_size);
                }                

                // clean frequency variables of outgoing edges in Wi 
                memset(Wi_freq, 0, sizeof(int)*255);   
                memset(DA_freq, 0, sizeof(int)*samples*255);
                memset(Wi_first_occurrence, 0, sizeof(int)*samples*255);
            }
            // if next LCP value is smaller than k-1 we have a new (k-1)-mer to keep track, so we clean W_freq values
            if(LCP[block_pos+1] < k-1){
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
            fwrite(reduced_LCP, sizeof(short), Wi_size, boss_reduced_LCP_file);

            // clean buffers
            memset(last, 0, sizeof(short)*10);   
            memset(W, 0, sizeof(char)*10);   
            memset(Wm, 0, sizeof(short)*10);   
            memset(colors, 0, sizeof(short)*10);   
            memset(coverage, 0, sizeof(int)*10);   
            memset(reduced_LCP, 0, sizeof(short)*10);   

            for(j = 0; j < 10; j++) coverage[j] = 1;

            Wi_size = 0; 
        }
        block_pos++;
        bi++;
    }

    // fix C values
    C[1] = C['$'];
    C[2] = C['A'] + C[1];
    C[3] = C['C'] + C[2];
    C[4] = C['G'] + C[3];
    C[5] = C['N'] + C[4];
    C[0] = 0;

    // Print BOSS construction info
    char alphabet[6] = {'$', 'A', 'C', 'G', 'N', 'T'};
    char boss_info[FILE_PATH];
    
    sprintf(boss_info, "results/%s-%s_k_%d.boss-info", file1, file2, k);
    
    FILE *boss_info_file = fopen(boss_info, "w");
                
    fprintf(boss_info_file, "Boss construction of %s and %s genomes merge:\n", file1, file2);
    fprintf(boss_info_file, "C array:\n");
    for(j = 0; j < 6; j++)
        fprintf(boss_info_file, "%c %d\n", alphabet[j], C[j]);
    fprintf(boss_info_file, "\n");

    fprintf(boss_info_file, "Boss length: %ld\n\n", i);

    fprintf(boss_info_file, "Total coverage: %ld\n\n", total_coverage);

    // Close used files
    fclose(boss_info_file);
    
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
    fclose(boss_reduced_LCP_file);

    // free BOSS construction needed variables
    free(LCP); free(BWT); free(DA); 
    
    // free BOSS construction variables
    free(last); free(W); free(Wm); free(colors); free(coverage); free(reduced_LCP);

    return i;
};