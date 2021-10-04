/**
   \file estringue.h
   \brief Functions for C-strings.
   \internal Guilherme P. Telles, 2014.
**/


#ifndef ESTRINGUE_H
#define ESTRINGUE_H

char* astrcat(char** dest, size_t* dest_size, char* src);

void chomp(char* str);
size_t rtrim(char* str);

int split(char* s, char* delimiters, char*** substr, int* n);

unsigned long hashvalue(char *str);

long* count_bytes(char* file);
int remap_bytes(char* filein, char* fileout, unsigned char* map);
int shift_alphabet(char* s, unsigned n, char ftag);

#endif
