#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "boss.h"

void Wi_sort(char *Wi, int *Wm, int *colors, int *coverage, int start, int end){
    int i;
    int range = end-start;
    char Wi_tmp[range];
    int Wm_tmp[range];
    int colors_tmp[range];
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

void add_edge(int i, char *W, int **last, int *colors, short *reduced_LCP, int freq, int *Wm, char bwt, int da, short lcp, int Wi_size, int edge_status){
    *W = bwt;
    *colors = da;
    *reduced_LCP = lcp;
    if(edge_status == 0){
        if(Wi_size == 0){
            (*last)[i] = 1;
        } else {
            (*last)[i-1] = 0;
            (*last)[i] = 1;
        }
    } else if(edge_status == 1){
        (*last)[i] = 1;
    } else if(edge_status == 2){
        (*last)[i-1] = 0;
        (*last)[i] = 1;
    }

    if(freq == 0){
        *Wm = 1;
    }
}

int boss_construction(FILE *mergeLCP, FILE *mergeDA, FILE *mergeBWT, int *C, int *last, char *W, int *Wm, int *colors, size_t n, int k, int samples, short *reduced_LCP, int *coverage, int *total_coverage, size_t docsSeparator, int mem){
    size_t i = 0; // iterates through Wi
    int j = 0;
    size_t bi = 0; // iterates through BWT, LCP and DA 
    int block_pos = 0; // iterates through BWT, LCP and DA memory blocks
    int Wi_size = 0; 
    int W_freq[255]; // frequency of outgoing edges in a (k-1)-mer suffix range (detects W- = 1)
    memset(W_freq, 0, sizeof(int)*255);
    int Wi_freq[255]; // frequency of outgoing edges in a k-mer suffix range (detects same outgoing edge in a vertex)
    memset(Wi_freq, 0,  sizeof(int)*255);
    int Wi_first_occurrence[samples][255]; // first occurence of an outgoing edge in a k-mer suffix range from a string collection
    memset(Wi_first_occurrence, 0,  sizeof(int)*samples*255);
    int DA_freq[samples][255]; // frequency of outgoing edges in a k-mer from a string collection (used to include same outgoing edge from distinct collections in BOSS representation)
    memset(DA_freq, 0,  sizeof(int)*samples*255);

    short *LCP = (short*)malloc((mem+1)*sizeof(short));
    int *DA = (int*)malloc(mem*sizeof(int));
    char *BWT = (char*)malloc(mem*sizeof(char));

    fread(LCP, sizeof(short), mem+1, mergeLCP);
    fread(DA, sizeof(int), mem, mergeDA);
    for(j = 0; j < mem; j++) DA[j] = DA[j] < docsSeparator ? 0 : 1;
    fread(BWT, sizeof(char), mem, mergeBWT);
    for(j = 0; j < mem; j++) BWT[j] = BWT[j] == 0 ? '$' : BWT[j];

    while(bi < n){
        // more than one outgoing edge of vertex i
        if(LCP[block_pos+1] >= k && bi != n-1){
            // since there is more than one outgoing edge, we don't need to check if BWT = $ or there is already BWT[bi] in Wi range
            if(BWT[block_pos] != '$'){ //change BWT[block_pos] with a new last_BWT variable
                if(Wi_freq[BWT[block_pos]] == 0){
                    // Add values to BOSS representation
                    add_edge(i, &W[i], &last, &colors[i], &reduced_LCP[i], W_freq[BWT[block_pos]], &Wm[i], BWT[block_pos], DA[block_pos], LCP[block_pos], Wi_size, 0);
                    Wi_first_occurrence[DA[block_pos]][BWT[block_pos]] = block_pos;
                    // Increment variables
                    C[BWT[block_pos]]++; W_freq[BWT[block_pos]]++; Wi_freq[BWT[block_pos]]++; DA_freq[DA[block_pos]][BWT[block_pos]]++; Wi_size++; i++;
                    (*total_coverage)++;
                } else {
                    // check if there is already outgoing edge labeled with BWT[bi] from DA[bi] leaving vertex i
                    if(DA_freq[DA[block_pos]][BWT[block_pos]] == 0){
                        add_edge(i, &W[i], &last, &colors[i], &reduced_LCP[i], W_freq[BWT[block_pos]], &Wm[i], BWT[block_pos], DA[block_pos], LCP[block_pos], Wi_size, 0);
                        Wi_first_occurrence[DA[block_pos]][BWT[block_pos]] = block_pos;
                        C[BWT[block_pos]]++; W_freq[BWT[block_pos]]++; Wi_freq[BWT[block_pos]]++; DA_freq[DA[block_pos]][BWT[block_pos]]++; Wi_size++; i++; 
                        (*total_coverage)++;
                    } else {
                        // increases the coverage information of the node with outgoing edge labeled with BWT[bi] from DA[bi] that is already on BOSS construction 
                        int existing_pos = Wi_first_occurrence[DA[block_pos]][BWT[block_pos]];
                        coverage[existing_pos]++;
                        (*total_coverage)++;
                    }
                }
            }
        } else {
            // just one outgoing edge of vertex i
            if(Wi_size == 0){
                add_edge(i, &W[i], &last, &colors[i], &reduced_LCP[i], W_freq[BWT[block_pos]], &Wm[i], BWT[block_pos], DA[block_pos], LCP[block_pos], Wi_size, 1);
                C[BWT[block_pos]]++; W_freq[BWT[block_pos]]++; i++;
                (*total_coverage)++;
            } 
            // last outgoing edge of vertex i
            else {
                // check if there is already outgoing edge labeled with BWT[bi] leaving vertex i
                if(Wi_freq[BWT[block_pos]] == 0){
                    add_edge(i, &W[i], &last, &colors[i], &reduced_LCP[i], W_freq[BWT[block_pos]], &Wm[i], BWT[block_pos], DA[block_pos], LCP[block_pos], Wi_size, 2);
                    C[BWT[block_pos]]++; W_freq[BWT[block_pos]]++; Wi_size++; i++;
                    (*total_coverage)++;
                } else {
                    // check if there is already outgoing edge labeled with BWT[bi] from DA[bi] leaving vertex i
                    if(DA_freq[DA[block_pos]][BWT[block_pos]] == 0){
                        add_edge(i, &W[i], &last, &colors[i], &reduced_LCP[i], W_freq[BWT[block_pos]], &Wm[i], BWT[block_pos], DA[block_pos], LCP[block_pos], Wi_size, 2);
                        C[BWT[block_pos]]++; W_freq[BWT[block_pos]]++; Wi_freq[BWT[block_pos]]++; DA_freq[DA[block_pos]][BWT[block_pos]]++; Wi_size++; i++;                     
                        (*total_coverage)++;
                    } else {
                        // increases the coverage information of the node with outgoing edge labeled with BWT[bi] from DA[bi] that is already on BOSS construction 
                        int existing_pos = Wi_first_occurrence[DA[block_pos]][BWT[block_pos]];
                        coverage[existing_pos]++;
                        (*total_coverage)++;
                    }
                }
                // sort outgoing edges of vertex i in lexigraphic order
                if(Wi_size > 1){
                    Wi_sort(W, Wm, colors, coverage, i-Wi_size, i);
                }
                // clean frequency variables of outgoing edges in Wi 
                memset(Wi_freq, 0, sizeof(int)*255);   
                memset(DA_freq, 0, sizeof(int)*samples*255);
                memset(Wi_first_occurrence, 0, sizeof(int)*samples*255);
            }
            // if next LCP value is smaller than k-1 we have a new (k-1)-mer to keep track, so we clean W_freq values
            if(LCP[bi+1] < k-1){
                memset(W_freq, 0,  sizeof(int)*255);
            }
            Wi_size = 0; 
        }

        if(bi%mem){
            // fseek(mergeLCP, bi, SEEK_SET);
            fread(LCP, sizeof(short), mem+1, mergeLCP);
            fread(DA, sizeof(int), mem, mergeDA);
            for(j = 0; j < mem; j++) DA[j] = DA[j] < docsSeparator ? 0 : 1;
            fread(BWT, sizeof(char), mem, mergeBWT);
            for(j = 0; j < mem; j++) BWT[j] = BWT[j] == 0 ? '$' : BWT[j];
            block_pos = 0;
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

void print_boss_result(int boss_len, char *file1, char *file2, int *C, int *last, char *W, int *Wm, int *colors, short *reduced_LCP, int *coverage, int total_coverage){
    int i;
    char alphabet[6] = {'$', 'A', 'C', 'G', 'N', 'T'};
    char boss_result[128];
    
    sprintf(boss_result, "results/%s-%s.boss", file1, file2);
    
    FILE *boss_file = fopen(boss_result, "w");
                
    fprintf(boss_file, "Boss construction of %s and %s genomes merge:\n", file1, file2);
    fprintf(boss_file, "C array:\n");
    for(i = 0; i < 6; i++)
        fprintf(boss_file, "%c %d\n", alphabet[i], C[i]);
    fprintf(boss_file, "\n");

    fprintf(boss_file, "Total coverage: %d\n\n", total_coverage);
    
    fprintf(boss_file, "BOSS:\nlast\tW\tW-\tcolor\tLCP\tcoverage\n");
    for(i = 0; i < boss_len; i++)
        fprintf(boss_file, "%4d\t%c\t%d\t%5d\t%3d\t%8d\n", last[i], W[i], Wm[i], colors[i], reduced_LCP[i], coverage[i]);

    fclose(boss_file);
}
