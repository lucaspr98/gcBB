#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "bwsd.h"

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
            double frac = (double)t[i]/(double)s;
            value += frac*log2(frac);
        }
    }

    return value*(-1);

}

void bwsd(int *DA, int n, double *expectation, double *entropy){
    int i;

    int *run_length = (int*)malloc((n*2)*sizeof(int));
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

    int *t = (int*)malloc(n*sizeof(int));
    memset(t, 0, n*sizeof(int));
    for(i = 0; i < pos; i+=2){
        t[run_length[i+1]]++;
    }

    *expectation = bwsd_expectation(t, pos/2, n);
    *entropy = bwsd_shannon_entropy(t, pos/2, n);
}

void print_bwsd_matrixes(double **Dm, double **De, char **files, int files_n){
    int i,j;
    FILE *bwsd_matrixes = fopen("results/bwsd_matrixes.txt", "w");
    fprintf(bwsd_matrixes, "Expectation matrix (D_m):\n");

    fprintf(bwsd_matrixes, "\t");
    for(i = 0; i < files_n; i++)
        fprintf(bwsd_matrixes, "%10d\t", i+1);
    fprintf(bwsd_matrixes, "\n");

    for(i = 0; i < files_n; i++){
        fprintf(bwsd_matrixes, "%d\t", i+1);
        for(j = 0; j < files_n; j++){
            if(j >= i)
                fprintf(bwsd_matrixes, "%10.5lf\t", Dm[i][j]);
            else
                fprintf(bwsd_matrixes, "%10.5lf\t", Dm[j][i]);
        }
        fprintf(bwsd_matrixes, "\n");
    }
    fprintf(bwsd_matrixes, "\n\n");

    fprintf(bwsd_matrixes, "Shannon's entropy matrix (D_e):\n");

    fprintf(bwsd_matrixes, "\t");
    for(i = 0; i < files_n; i++)
        fprintf(bwsd_matrixes, "%10d\t", i+1);
    fprintf(bwsd_matrixes, "\n");

    for(i = 0; i < files_n; i++){
        fprintf(bwsd_matrixes, "%d\t", i+1);
        for(j = 0; j < files_n; j++){
            if(j >= i)
                fprintf(bwsd_matrixes, "%10.5lf\t", De[i][j]);
            else
                fprintf(bwsd_matrixes, "%10.5lf\t", De[j][i]);
        }
        fprintf(bwsd_matrixes, "\n");
    }
    
    fprintf(bwsd_matrixes, "\n\n");
    fprintf(bwsd_matrixes, "Matrix info:\n");
    for(i = 0; i < files_n; i++)
        fprintf(bwsd_matrixes, "%d: %s\n", i+1, files[i]);

    fprintf(bwsd_matrixes, "\n\n");
    fprintf(bwsd_matrixes, "CSV expectation distance matrix:\n");

    for(i = 0; i < files_n; i++)
        fprintf(bwsd_matrixes, ",%s",files[i]);
    fprintf(bwsd_matrixes, "\n");
    for(i = 0; i < files_n; i++){
        fprintf(bwsd_matrixes, "%s,", files[i]);
        for(j = 0; j < files_n; j++){
            if(j >= i)
                fprintf(bwsd_matrixes, "%lf,", Dm[i][j]);
            else
                fprintf(bwsd_matrixes, "%lf,", Dm[j][i]);
        }
        fprintf(bwsd_matrixes, "\n");
    }
    
    fprintf(bwsd_matrixes, "\n\n");
    fprintf(bwsd_matrixes, "CSV entropy distance matrix:\n");

    for(i = 0; i < files_n; i++)
        fprintf(bwsd_matrixes, ",%s",files[i]);
    fprintf(bwsd_matrixes, "\n");
    for(i = 0; i < files_n; i++){
        fprintf(bwsd_matrixes, "%s,", files[i]);
        for(j = 0; j < files_n; j++){
            if(j >= i)
                fprintf(bwsd_matrixes, "%lf,", De[i][j]);
            else
                fprintf(bwsd_matrixes, "%lf,", De[j][i]);
        }
        fprintf(bwsd_matrixes, "\n");
    }

    fprintf(bwsd_matrixes, "\n\n");
    fprintf(bwsd_matrixes, "Phylip expectation distance matrix:\n");

    fprintf(bwsd_matrixes, "%d\n", files_n);
    for(i = 0; i < files_n; i++){
        fprintf(bwsd_matrixes, "%s\t", files[i]);
        for(j = 0; j < files_n; j++){
            if(j >= i)
                fprintf(bwsd_matrixes, "\t%lf", Dm[i][j]);
            else
                fprintf(bwsd_matrixes, "\t%lf", Dm[j][i]);
        }
        fprintf(bwsd_matrixes, "\n");
    }

    fprintf(bwsd_matrixes, "\n\n");
    fprintf(bwsd_matrixes, "Phylip entropy distance matrix:\n");

    fprintf(bwsd_matrixes, "%d\n", files_n);
    for(i = 0; i < files_n; i++){
        fprintf(bwsd_matrixes, "%s\t", files[i]);
        for(j = 0; j < files_n; j++){
            if(j >= i)
                fprintf(bwsd_matrixes, "\t%lf", De[i][j]);
            else
                fprintf(bwsd_matrixes, "\t%lf", De[j][i]);
        }
        fprintf(bwsd_matrixes, "\n");
    }

    fclose(bwsd_matrixes);
}