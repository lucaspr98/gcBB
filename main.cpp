#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void boss_construction(int *LCP, int *DA, char *BWT, int *C, int *last, int *Wm, unsigned char *W, int n, int k);

void Wi_sort(unsigned char *Wi, int *Wm, int start, int end);

void boss_example();

void colored_boss_example();

int main(){

    // boss_example();

    // colored_boss_example();

}

void boss_example() {
    int n0, n1, k;
    n0 = 14;
    n1 = 7;
    k = 3;

    char BWT_0[n0+1] = "TTC$CCTATAAA$C";
    int LCP_0[n0] = {-1,3,0,2,1,3,0,2,3,1,0,3,1,3};

    char BWT_1[n1+1] = "GCTGA$C";
    int LCP_1[n1] = {-1,0,0,1,0,1,0};
    
    /*** Construct BOSS representation ***/

    int *last_0 = (int*)malloc((n0)*sizeof(int*));
    int *Wm_0 = (int*)malloc((n0)*sizeof(int*));

    int *last_1 = (int*)malloc((n1)*sizeof(int*));
    int *Wm_1 = (int*)malloc((n1)*sizeof(int*));

    // create and initialize C map
    int C_0[255];
    memset(C_0, 0, sizeof(int)*255);
    int C_1[255];
    memset(C_1, 0, sizeof(int)*255);

    // create and initizalize W string
    unsigned char *W_0 = (unsigned char*)malloc((n0)*sizeof(unsigned char*));
    unsigned char *W_1 = (unsigned char*)malloc((n1)*sizeof(unsigned char*));

    // BOSS construction
    boss_construction(LCP_0, NULL, BWT_0, C_0, last_0, Wm_0, W_0, n0, k);
    boss_construction(LCP_1, NULL, BWT_1, C_1, last_1, Wm_1, W_1, n1, k);
}

void colored_boss_example() {
    int n, k;
    n = 21;
    k = 3;

    char BWT[n+1] = "TTGC$CCCTTATGAA$AA$CC";
    int LCP[n] = {-1,3,3,0,2,1,1,3,0,2,2,3,1,4,0,1,0,3,1,3,3};
    int DA[n] = {0,0,1,0,0,1,0,0,0,1,0,0,1,0,1,1,0,0,0,1,0};
    
    /*** Construct BOSS representation ***/

    int *last = (int*)malloc((n)*sizeof(int*));
    int *Wm = (int*)malloc((n)*sizeof(int*));

    // create and initialize C map
    int C[255];
    memset(C, 0, sizeof(int)*255);

    // create and initizalize W string
    unsigned char *W = (unsigned char*)malloc((n)*sizeof(unsigned char*));

    // BOSS construction
    boss_construction(LCP, DA, BWT, C, last, Wm, W, n, k);
}

void Wi_sort(unsigned char *Wi, int *Wm, int start, int end){
    int i;
    int range = end-start;
    int Wi_tmp[range];
    int Wi_aux[255];
    int Wm_tmp[range];
    int Wm_aux[2] = {0,0};

    memset(Wi_tmp, 0, sizeof(int)*(range));    
    memset(Wm_tmp, 0, sizeof(int)*(range)); 
    memset(Wi_aux, 0, sizeof(int)*255);

    for(i = start; i < end; i++){
        Wi_aux[Wi[i]]++;
    }

    for(i = 1; i < 255; i++){
        Wi_aux[i] += Wi_aux[i-1];
    }

    Wm_aux[1] += Wm_aux[0];

    for(i = start; i < end; i++){
        Wi_tmp[Wi_aux[Wi[i]-1]] = Wi[i];
        Wm_tmp[Wi_aux[Wi[i]-1]] = Wm[i];
        Wi_aux[Wi[i]]--;
    }

    for(i = start; i < end; i++){
        Wi[i] = Wi_tmp[i-start];
        Wm[i] = Wm_tmp[i-start];
    }
}

void boss_construction(int *LCP, int *DA, char *BWT, int *C, int *last, int *Wm, unsigned char *W, int n, int k){
    int i = 0; // iterates through Wi
    int j = 0; // auxiliary iterator  
    int bi = 0; // iterates through BWT and LCP
    int Wi_size = 0; 
    int W_freq[255]; // frequency of outgoing edges in a (k-1)-mer suffix range (detect W- = 1)
    memset(W_freq, 0, sizeof(int)*255);
    int Wi_freq[255]; // frequency of outgoing edges in a k-mer suffix range (detect same outgoing edge in a vertex)
    memset(Wi_freq, 0,  sizeof(int)*255);
    int DA_freq[255]; // frequency of outgoing edges in a k-mer from a string collection (include same outgoing edge from distinct collections)
    memset(Wi_freq, 0,  sizeof(int)*255);

    while(bi < n-1){
        // more than one outgoing edge of vertex i
        if(LCP[bi+1] >= k){
            // since there is more than one outgoing edge, we don't need to check if BWT = $ or there is already BWT[bi] in Wi range
            if(BWT[bi] != '$' && Wi_freq[BWT[bi]] == 0){
                W[i] = BWT[bi];
                if(Wi_size == 0){
                    last[i] = 1;
                } else {
                    last[i-1] = 0;
                    last[i] = 1;
                }
                if(W_freq[BWT[bi]] == 0){
                    Wm[i] = 1;
                }
                C[BWT[bi]]++;
                W_freq[BWT[bi]]++;
                Wi_freq[BWT[bi]]++;
                Wi_size++;
                i++; 
            }
        } else {
            // just one outgoing edge of vertex i
            if(Wi_size == 0){
                W[i] = BWT[bi];
                last[i] = 1;
                if(W_freq[BWT[bi]] == 0){
                    Wm[i] = 1;
                }
                W_freq[BWT[bi]]++;
                i++;
                C[BWT[bi]]++;
            } 
            // last outgoing edge of vertex i
            else {
                // check if there is already outgoing edge labeled with BWT[bi] leaving vertex i
                if(Wi_freq[BWT[bi]] == 0){
                    W[i] = BWT[bi];
                    last[i-1] = 0;
                    last[i] = 1;
                    if(W_freq[BWT[bi]] == 0){
                        Wm[i] = 1;
                    }
                    i++;
                    Wi_size++;
                    C[BWT[bi]]++;
                    W_freq[BWT[bi]]++;
                } 
                // sort outgoing edges of vertex i in lexigraphic order
                if(Wi_size > 1){
                    Wi_sort(W, Wm, i-Wi_size, i);
                }
                memset(Wi_freq, 0, sizeof(int)*255);   
            }
            // if next LCP value is smaller than k-1 we have a new (k-1)-mer to keep track, so we clean W_freq values
            if(LCP[bi+1] < k-1){
                memset(W_freq, 0,  sizeof(int)*255);
            }
            Wi_size = 0; 
        }
        
        bi++;
    }
    // check last element
    if(Wi_size == 0){
        W[i] = BWT[bi];
        last[i] = 1;
        if(W_freq[BWT[bi]] == 0){
            Wm[i] = 1;
        }
        C[BWT[bi]]++;
        i++;
    } else {
        if(Wi_freq[BWT[bi]] == 0){
            W[i] = BWT[bi];
            last[i-1] = 0;
            last[i] = 1;
            if(W_freq[BWT[bi]] == 0){
                Wm[i] = 1;
            }
            i++;
            Wi_size++;
            C[BWT[bi]]++;
            W_freq[BWT[bi]]++;
        }
        if(Wi_size > 1) 
            Wi_sort(W, Wm, i-Wi_size, i);
    }

    // fix C values
    C[1] = C['$'];
    C[2] = C['A'] + C[1];
    C[3] = C['C'] + C[2];
    C[4] = C['G'] + C[3];

    C[0] = 0;

    char alphabet[5] = {'$', 'A', 'C', 'G', 'T'};
    // Print results 
    printf("C array:\n");
    for(j = 0; j < 5; j++){
        printf("%c %d\n", alphabet[j], C[j]);
    }
    printf("\n");

    printf("BOSS:\nlast\tW\tW-\n");
    for(j = 0; j < i; j++){
        printf("%d\t%c\t%d\n", last[j], W[j], Wm[j]);
    }
    printf("\n");
};
