#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "boss.h"

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
