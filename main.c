#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void compute_files(char *file1, char *file2, int k);

void boss_construction(int *LCP, int *DA, char *BWT, int *C, int *last, char *W, int *Wm, int *colors, int n, int k, int samples);

void Wi_sort(char *Wi, int *Wm,  int *colors, int start, int end);

void bwsd(int *DA, int n, double *expectation, double *entropy);

double bwsd_expectation(int *t, int s, int n);

double bwsd_shannon_entropy(int *t, int s, int n);

// void boss_example(); // do we want to compute boss to a single genome? If we want we need a few ifdef

void colored_boss_example();

int main(int argc, char *argv[]){
    int i;

    // colored_boss_example();
    // exit(-1);

    /******** Check arguments and compute external needed files ********/

    // Check arguments
    if (argc < 4) {
        printf("Missing argument!\nusage: ./bgc [file1.fastq] [file2.fastq] [k]\n");
        exit(-1);
    }

    char *file1 = argv[1];
    char *file2 = argv[2];
    int k = atoi(argv[3]);
    
    // // Computes BWT, LCP and DA from both files
    compute_files(file1, file2, k);

    FILE *mergeBWT = fopen("merge.bwt", "r");
    FILE *mergeLCP = fopen("merge.4.lcp", "rb");
    FILE *mergeDA = fopen("merge.4.da", "rb");

    char docs[128];
    sprintf(docs, "%s.docs", file1);
    FILE *mergeDocs = fopen(docs, "rb");    

    fseek(mergeBWT, 0, SEEK_END);
    int n = ftell(mergeBWT);
    rewind(mergeBWT);

    int docsSize;
    fread(&docsSize, 4, 1, mergeDocs);

    // /******** Construct BOSS representation ********/

    // Declare needed variables
    char *BWT;
    int *LCP, *DA;
    int *last, *Wm, *colors;
    char *W;
    int C[255];
    int samples = 2;

    // Initialize variables needed to construct BOSS
    BWT = (char*)malloc(n*sizeof(char));
    LCP = (int*)malloc(n*sizeof(int));
    DA = (int*)malloc(n*sizeof(int));

    fread(BWT, sizeof(char), n, mergeBWT);
    fread(LCP, sizeof(int), n, mergeLCP);
    int buff;
    for(i = 0; i < n; i++){
        fread(&buff, 4, 1, mergeDA);
        if(buff < docsSize)
            DA[i] = 0; 
        else
            DA[i] = 1;

        if(BWT[i] == 0)
            BWT[i] = '$';
    }

    // Initialize BOSS variables
    memset(C, 0, sizeof(int)*255);
    last = (int*)malloc(n*sizeof(int));
    Wm = (int*)malloc(n*sizeof(int));
    colors = (int*)malloc(n*sizeof(int));
    W = (char*)malloc(n*sizeof(char));

    boss_construction(LCP, DA, BWT, C, last, W, Wm, colors, n, k, samples);

    // BWSD
    double expectation, entropy;

    bwsd(DA, n, &expectation, &entropy);

}

void compute_files(char *file1, char *file2, int k){
    char eGap1[128];
    char eGap2[128];
    char eGapMerge[128];

    sprintf(eGap1, "egap/eGap --trlcp %d --lbytes 4 --da %s", k, file1);
    sprintf(eGap2, "egap/eGap --trlcp %d --lbytes 4 --da %s", k, file2);
    sprintf(eGapMerge, "egap/eGap --bwt --trlcp %d --lbytes 4 --da -o merge %s.bwt %s.bwt", k, file1, file2);

    system(eGap1);
    system(eGap2);
    system(eGapMerge);
}

// void boss_example() {
//     int n0, n1, k;
//     n0 = 14;
//     n1 = 7;
//     k = 3;

//     char BWT_0[n0+1];
//     strcpy(BWT_0, "TTC$CCTATAAA$C");
//     int LCP_0[] = {-1,3,0,2,1,3,0,2,3,1,0,3,1,3};

//     char BWT_1[n1+1]; 
//     strcpy(BWT_1, "GCTGA$C");
//     int LCP_1[] = {-1,0,0,1,0,1,0};
    
//     /*** Construct BOSS representation ***/

//     int *last_0 = (int*)malloc((n0)*sizeof(int*));
//     int *Wm_0 = (int*)malloc((n0)*sizeof(int*));

//     int *last_1 = (int*)malloc((n1)*sizeof(int*));
//     int *Wm_1 = (int*)malloc((n1)*sizeof(int*));

//     // create and initialize C map
//     int C_0[255];
//     memset(C_0, 0, sizeof(int)*255);
//     int C_1[255];
//     memset(C_1, 0, sizeof(int)*255);

//     // create and initizalize W string
//     char *W_0 = (char*)malloc((n0)*sizeof(char*));
//     char *W_1 = (char*)malloc((n1)*sizeof(char*));

//     // BOSS construction
//     boss_construction(LCP_0, NULL, BWT_0, C_0, last_0, W_0, Wm_0, NULL, n0, k, 1);
//     boss_construction(LCP_1, NULL, BWT_1, C_1, last_1, W_1, Wm_1, NULL, n1, k, 1);
// }

void colored_boss_example() {
    int n, k;
    n = 23;
    k = 3;

    char BWT[n+1];
    strcpy(BWT, "TTGC$CCCTTATGAA$AA$CCTT");
    int LCP[] = {-1,3,3,0,2,1,1,3,0,2,2,3,1,4,0,1,0,3,1,3,3,4,6};
    int DA[] = {0,0,1,0,0,1,0,0,0,1,0,0,1,0,1,1,0,0,0,1,0,1,0};
    
    /*** Construct BOSS representation ***/

    int *last = (int*)malloc((n)*sizeof(int*));
    int *Wm = (int*)malloc((n)*sizeof(int*));
    int *colors = (int*)malloc((n)*sizeof(int*));

    // create and initialize C map
    int C[255];
    memset(C, 0, sizeof(int)*255);

    // create and initizalize W string
    char *W = (char*)malloc((n)*sizeof(char*));

    // BOSS construction
    boss_construction(LCP, DA, BWT, C, last, W, Wm, colors, n, k, 2);
}

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

void boss_construction(int *LCP, int *DA, char *BWT, int *C, int *last, char *W, int *Wm, int *colors, int n, int k, int samples){
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
            // if(BWT[bi] == 0)
            //     BWT[bi] = '$';
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

    char alphabet[6] = {'$', 'A', 'C', 'G', 'N', 'T'};
    // Print results 
    printf("C array:\n");
    for(j = 0; j < 6; j++){
        printf("%c %d\n", alphabet[j], C[j]);
    }
    printf("\n");

    printf("BOSS:\nlast\tW\tW-\tcolor\n");
    for(j = 0; j < i; j++){
        printf("%d\t\t%c\t%d\t%d\n", last[j], W[j], Wm[j], colors[j]);
    }
};

double bwsd_expectation(int *t, int s, int n){
    int i;
	double value = 0.0;
    double frac;

	
	for(i = 1; i < n; i++){
        if(t[i] != 0){
            frac = (double)t[i]/s;
            value += i*frac;
        }
    }

    return value-1;
}

double bwsd_shannon_entropy(int *t, int s, int n){
    int i;
	double value = 0.0;
	
    for(i = 1; i < n; i++){
        if(t[i] != 0){
            double frac = (double)t[i]/s;
            value += frac*log2(frac);
        }
    }

    return value*(-1);

}

void bwsd(int *DA, int n, double *expectation, double *entropy){
    int i;

    int *run_length = (int*)malloc((n+6)*sizeof(int));
    int current = 0;
    run_length[0] = 0;
    run_length[1] = 0;
    int pos = 1;
    for(i = 0; i < n; i++){
        if(DA[i] == current)
            run_length[pos]++;
        else {
            current = DA[i];
            run_length[pos+1]=current;
            run_length[pos+2]=1;
            pos += 2;
        }
    }
    if(run_length[pos-1] == 0){
        pos+= 2;
        run_length[pos-1] = 1;
        run_length[pos] = 0;
    }
    pos++; //size of run_lentgh

    int t[n];
    memset(t, 0, (n)*sizeof(int));
    for(i = 0; i < pos; i+=2){
        t[run_length[i+1]]++;
    }

    *expectation = bwsd_expectation(t, pos/2, n);
    *entropy = bwsd_shannon_entropy(t, pos/2, n);

    printf("BWSD:\n");
    printf("Expectation: %.5lf\nShannons Entropy: %.5lf\n", *expectation, *entropy);
}

