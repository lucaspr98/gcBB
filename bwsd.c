#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <libgen.h>
#include "bwsd.h"

#ifndef COVERAGE
	#define COVERAGE 0
#endif

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

void bwsd(int *colors, short *reduced_LCP, int *coverage, int n, int k, double *expectation, double *entropy, int mem){
    int i;

    int *run_length = (int*)malloc((n*3)*sizeof(int));
    int current = 0;
    run_length[0] = 0;
    run_length[1] = 0;
    int pos = 1;

    for(i = 0; i < n; i++){
        #if COVERAGE 
        if(reduced_LCP[i] > k && reduced_LCP[i+1] > k && colors[i] != colors[i+1]){
            if(coverage[i] >= coverage[i+1]){
                int repetitions = 0;
                int flip = colors[i];
                while(repetitions <= coverage[i+1]){
                    run_length[pos+1] = flip;
                    run_length[pos+2] = 1;
                    pos += 2;
                    flip = flip == colors[i] ? colors[i+1] : colors[i];
                    repetitions++;
                }
                int remaining = coverage[i] - repetitions;
                run_length[pos+1] = colors[i];
                run_length[pos+2] = remaining;
                pos += 2; 
            } 
            else {
                int repetitions = 0;
                int flip = colors[i];
                while(repetitions <= coverage[i]){
                    run_length[pos+1] = flip;
                    run_length[pos+2] = 1;
                    pos += 2;
                    flip = flip == colors[i] ? colors[i+1] : colors[i];
                    repetitions++;
                }
                int remaining = coverage[i+1] - repetitions;
                run_length[pos+1] = colors[i+1];
                run_length[pos+2] = remaining;
                pos += 2; 
            }
        } else { 
        #endif
            if(colors[i] == current)
                run_length[pos]++;
            else {
                current = colors[i];
                run_length[pos+1]=current;
                run_length[pos+2]=1;
                pos += 2;
            }
        #if COVERAGE
        }
        #endif

        i++;
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

void print_bwsd_matrixes(double **Dm, double **De, char **files, int files_n, char *path){
    int i,j;
    char *ptr;
    char outputFile[128];

    int len = strlen(path);
    char folder[len];
    strcpy(folder, basename(path));

    if(files_n > 2){
        ptr = strchr(path, '/');
        if (ptr != NULL)
            *ptr = '\0';

        #if COVERAGE
            sprintf(outputFile, "results/%s_distance_matrixes_coverage_1.txt", folder);
        #else
            sprintf(outputFile, "results/%s_distance_matrixes_coverage_0.txt", folder);
        #endif
    } else {
        #if COVERAGE
            sprintf(outputFile, "results/%s-%s_distance_matrixes_coverage_1.txt", files[0], files[1]);
        #else
            sprintf(outputFile, "results/%s-%s_distance_matrixes_coverage_0.txt", files[0], files[1]);
        #endif
    }
    

    FILE *bwsd_matrixes = fopen(outputFile, "w");
    
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