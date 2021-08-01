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

size_t apply_coverage(short primaryColor, short secondaryColor, int primaryCoverage, int secondaryCoverage, short *rl_color, int *rl_freq, size_t pos){
    int repetitions = 0;
    int flip = primaryColor;
    while(repetitions <= secondaryCoverage){
        rl_color[pos] = flip;
        rl_freq[pos] = 1;
        pos++;
        flip = flip == primaryColor ? secondaryColor : primaryColor;
        repetitions++;
    }
    int remaining = primaryCoverage - repetitions;
    // rl_color[pos] = primaryColor;
    // rl_freq[pos] = remaining;

    return pos++;
}

void bwsd(char* file1, char* file2, size_t n, int k, double *expectation, double *entropy, int mem){
    size_t i;

    size_t size = n+1;
    
    short *rl_color = (short*)malloc(size*sizeof(short));
    int *rl_freq = (int*)malloc(size*sizeof(int));

    int current = 0;
    rl_color[0] = 0;
    rl_freq[0] = 0;
    size_t pos = 0;

    short *colors = (short*)malloc((mem+2)*sizeof(short));
    short *reduced_LCP = (short*)malloc((mem+2)*sizeof(short));
    int *coverage = (int*)malloc((mem+2)*sizeof(int));

    char color_file_name[128];
    char reduced_LCP_file_name[128];
    char coverage_file_name[128];

    sprintf(color_file_name, "results/%s-%s.2.colors", file1, file2);
    sprintf(reduced_LCP_file_name, "results/%s-%s.2.reduced_lcp", file1, file2);
    sprintf(coverage_file_name, "results/%s-%s.4.coverage", file1, file2);
    
    FILE *colors_file = fopen(color_file_name, "rb");
    FILE *recuced_LCP_file = fopen(reduced_LCP_file_name, "rb");
    FILE *coverage_file = fopen(coverage_file_name, "rb");

    fread(colors, sizeof(short), mem+1, colors_file);
    fread(reduced_LCP, sizeof(short), mem+1, colors_file);
    fread(coverage, sizeof(int), mem+1, coverage_file);

    int block_pos = 0;

    for(i = 0; i < n; i++){

        if(i != 0 && i%mem == 0){
            fseek(colors_file, -2L, SEEK_CUR);
            fread(colors, sizeof(short), mem+1, colors_file);

            fseek(recuced_LCP_file, -2L, SEEK_CUR);
            fread(reduced_LCP, sizeof(short), mem+1, recuced_LCP_file);

            fseek(coverage_file, -4L, SEEK_CUR);
            fread(coverage, sizeof(int), mem+1, coverage_file);

            block_pos=0;
        }

        #if COVERAGE 
        if(reduced_LCP[block_pos] > k && reduced_LCP[block_pos+1] > k && colors[block_pos] != colors[block_pos+1] && coverage[block_pos] > 1 ){
            if(coverage[block_pos] >= coverage[block_pos+1]){
                pos = apply_coverage(colors[block_pos], colors[block_pos+1], coverage[block_pos], coverage[block_pos+1], rl_color, rl_freq, pos);
            } 
            else {
                pos = apply_coverage(colors[block_pos], colors[block_pos+1], coverage[block_pos+1], coverage[block_pos], rl_color, rl_freq, pos);
            }
            pos++;
        } else { 
        #endif
            if(colors[block_pos] == current){
                rl_freq[pos]++;
            } else {
                current = colors[block_pos];
                pos++;
                rl_color[pos]=current;
                rl_freq[pos]=1;
            }
        #if COVERAGE
        }
        #endif

        block_pos++;
    }

    free(colors); free(coverage); free(reduced_LCP);

    if(rl_freq[pos] == 0){
        pos++;
        rl_color[pos] = 1;
        rl_freq[pos] = 0;
    }
    pos++; //size of run_lentgh

    int *t = (int*) calloc((n+1), sizeof(int));
    for(i = 0; i < pos; i++)
        t[rl_freq[i]]++;
    

    size_t s = pos/2;
    *expectation = bwsd_expectation(t, s, n);
    *entropy = bwsd_shannon_entropy(t, s, n);

    free(rl_color); free(rl_freq); free(t);

    fclose(colors_file);
    fclose(recuced_LCP_file);
    fclose(coverage_file);
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