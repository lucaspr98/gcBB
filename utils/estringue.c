/*
  Guilherme P. Telles.
*/

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>

#include "estringue.h"



/**
   \brief Trim the right end of a string.

   Remove the rightmost blanks (spaces, tabs, newlines) from a string.
**/
size_t rtrim(char* str) {
  size_t l = strlen(str)-1;
  while (l >= 0 && (str[l] == ' ' || str[l] == '\n' || str[l] == '\t' ||
                    str[l] == '\r' || str[l] == '\f' || str[l] == '\v'))
    l--;
  str[l+1] = 0;

  return l+1;
}




/**
   \brief Evaluate a hash value for a string.

   This is the function djb2.
**/
unsigned long hashvalue(char *str) {

  unsigned long hash = 5381;
  int c;

  while ((c = *str++))
    hash = ((hash << 5) + hash) + c; // hash * 33 + c

  return hash;
}




/**
   \brief Remove a single \n from a string.

   If the string terminates with \\n then \\n is replaced by \\0.
**/
void chomp(char* s) {

  int n = strlen(s);
  if (s[n-1] == '\n')
    s[n-1] = 0;
}




/**
   \brief Split a string at a delimiter.

   Split s at every occurrence of a delimiter.  Returns the address of a new
   array of substrings in *substr and the number of substrings in *n.  If no
   delimiter occurrs in s, then it sets *n to 1 and places a copy of s at
   substr[0].  Adding \0 as a delimiter will probably not work.

   \return On success it returns the number of substrings.  On failure it
   returns 0, sets *n to 0 and leaves errno as set by malloc() on failure.
**/
int split(char* s, char* delimiters, char*** substr, int* n) {

  int i,j,k;

  i = k = 0;
  while (s[i]) {
    if (index(delimiters,s[i]) != NULL)
      k++;
    i++;
  }

  *n = k+1;
  char** S = calloc(*n,sizeof(char*));
  if (!S) goto nomemh;

  if (k == 0) {
    S[0] = strdup(s);
    *substr = S;
    return 1;
  }


  i = k = 0;
  while (k<(*n)) {

    j = i;
    while (index(delimiters,s[j]) == NULL && s[j])
      j++;

    if (i==j) {
      S[k] = calloc(1,sizeof(char));
      if (!S[k]) goto nomemh;
    }
    else {
      S[k] = malloc((j-i+1)*sizeof(char));
      if (!S[k]) goto nomemh;
      memcpy(S[k],s+i,j-i*sizeof(char));
      S[k][j-i] = 0;
    }

    k++;
    i = j+1;
  }

  *substr = S;
  return (*n);

 nomemh:
  k = errno;
  while (--(*n) >= 0)
    free(S[*n]);
  free(S);
  *n = 0;
  errno = k;
  return 0;
}



/**
   \brief strcat with memory resizing.

   This function performs size checking and resizing of *dest as needed, and
   then copies src to the end of *dest.

   \return On success it returns the address of dest.  On failure it returns
   NULL and errno remains as set by realloc().
**/
char* astrcat(char** dest, size_t* dest_size, char* src) {

  size_t n;
  if (*dest)
    n = strlen(*dest);
  else
    n = 0;

  size_t m = strlen(src);

  if (*dest_size < n+m+1) {
    char* s = realloc(*dest,(n+m+1)*sizeof(char));
    if (!s) return 0;

    *dest = s;
    *dest_size = n+m+1;
  }

  return strcpy(*dest+n,src);
}



/**
   \brief Count the frequencies of each byte in a file.

   \return A new array with 256 counters.  On failure it returns 0 and errno
   remains as set by fopen() or is set to EIO if any error happens when the file
   is read.
**/
long* count_bytes(char* file) {

  long* count = calloc(256,sizeof(long));
  if (!count) return 0;

  FILE* f = fopen(file,"r");
  if (!f) { free(count); return 0; }

  int c;
  while ((c = fgetc(f)) != EOF)
    count[c]++;

  if (ferror(f)) {
    fclose(f);
    errno = EIO;
    free(count);
    return 0;
  }

  fclose(f);

  return count;
}



/**
   \brief Remap bytes in a file.

   Reads filein and creates a new file with bytes replaced according to map.  If
   the byte at position i of filein is b, then the byte at position i of fileout
   will be map[b].

   \return On success it returns 1.  On failure it returns 0 and errno remains
   as set by fopen() or is set to EIO if any error happens when the files are read
   or written.
**/
int remap_bytes(char* filein, char* fileout, unsigned char* map) {

  FILE* in = fopen(filein,"r");
  if (!in)  return 0;

  FILE* out = fopen(fileout,"w");
  if (!out)  { fclose(in); return 0; }

  int c;
  while ((c = fgetc(in)) != EOF && fputc(map[c],out) != EOF) ;

  if (ferror(in) || ferror(out)) {
    fclose(in);
    fclose(out);
    errno = EIO;
    return 0;
  }

  fclose(in);
  fclose(out);

  return 1;
}



/**
   \brief Shift the alphabet of a string.

   This function shifts the alphabet of a string such that the smaller char in s
   gets value ftag, the second smaller char gets value ftag+1 and so on.  If the
   input string has 128 distinct characters then this function does nothing.  If
   a char would overflow with the result of ftag+k then this function does
   nothing.

   \return On success it returns the size of the shifted alphabet.  On failure
   it returns 0 and leaves errno as set by malloc() or sets it to EOVERFLOW if
   an overflow would occur.  If the input string has 128 distinct characters
   then it returns 128.
**/
int shift_alphabet(char* s, unsigned n, char ftag) {

  unsigned char* mark = calloc(128,sizeof(unsigned char));
  if (!mark) return 0;

  int i;
  for (i=0; i<n; i++)
    if (!mark[s[i]])
      mark[s[i]] = 1;

  for (i=1; i<128; i++)
    mark[i] += mark[i-1];

  if (mark[127] == 128)
    return 128;

  if (((int)mark[127])+((int)ftag) > 127) {
    errno = EOVERFLOW;
    return 0;
  }

  for (i=0; i<n; i++)
    s[i] = mark[s[i]]-1+ftag;

  return mark[127];
}
