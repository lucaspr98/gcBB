#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct bucket_entry
{
    char *str;
    size_t len;
    struct bucket_entry *next;
};
typedef struct bucket_entry bucket_entry;
 
/* A linked list */
struct bucket
{
    bucket_entry *head;
    bucket_entry *tail;
};
typedef struct bucket bucket;
 
/* Initialise buckets */
static void init_buckets(bucket *buckets)
{
    unsigned int b;
    for (b = 0; b < 256; b++) {
        buckets[b].head = NULL;
        buckets[b].tail = NULL;
    }
}
 
/*
 * Turn entries into a linked list of words with their lengths
 * Return the length of the longest word
 */
static size_t init_entries(char **strings, size_t len, bucket_entry *entries)
{
    unsigned int s;
    size_t maxlen = 0;
    for (s = 0; s < len; s++) {
        entries[s].str = strings[s];
        entries[s].len = strlen(strings[s]);
        if (entries[s].len > maxlen) {
            maxlen = entries[s].len;
        }
        if (s < len - 1) {
            entries[s].next = &entries[s + 1];
        }
        else {
            entries[s].next = NULL;
        }
    }
    return maxlen;
}
 
/* Put strings into buckets according to the character c from the right */
void bucket_strings(bucket_entry *head, bucket *buckets, unsigned int c)
{
    bucket_entry *current = head;
    while (current != NULL) {
        bucket_entry * next = current->next;
        current->next = NULL;
        int pos = current->len - 1 - c;
        unsigned char ch;
        if (pos < 0) {
            ch = 0;
        }
        else {
            ch = current->str[pos];
        }
        if (buckets[ch].head == NULL) {
            buckets[ch].head = current;
            buckets[ch].tail = current;
        }
        else {
            buckets[ch].tail->next = current;
            buckets[ch].tail = buckets[ch].tail->next;
        }
        current = next;
    }
}
 
/* Merge buckets into a single linked list */
bucket_entry *merge_buckets(bucket *buckets)
{
    bucket_entry *head = NULL;
    bucket_entry *tail;
    unsigned int b;
    for (b = 0; b < 256; b++) {
        if (buckets[b].head != NULL) {
            if (head == NULL) {
                head = buckets[b].head;
                tail = buckets[b].tail;
            }
            else {
                tail->next = buckets[b].head;
                tail = buckets[b].tail;
            }
        }
    }
    return head;
}
 
void print_buckets(const bucket *buckets)
{
    unsigned int b;
    for (b = 0; b < 256; b++) {
        if (buckets[b].head != NULL) {
            const bucket_entry *entry;
            printf("[%d] ", b);
            for (entry = buckets[b].head; entry != NULL; entry = entry->next) {
                printf("%s", entry->str);
                if (entry->next) {
                    printf(" -> ");
                }
            }
            putchar('\n');
        }
    }
    putchar('\n');
}
 
void print_list(const bucket_entry *head)
{
    const bucket_entry *current;
    for (current = head; current != NULL; current = current->next) {
        printf("%s", current->str);
        if (current->next != NULL) {
            printf(" -> ");
        }
    }
    printf("\n");
}
 
void copy_list(const bucket_entry *head, char **strings)
{
    const bucket_entry *current;
    unsigned int s;
    for (current = head, s = 0; current != NULL; current = current->next, s++) {
        strings[s] = current->str;
    }
}
 
void radix_sort_string(char **strings, size_t len)
{
    size_t maxlen;
    unsigned int c;
    bucket_entry *head;
    bucket_entry *entries = malloc(len * sizeof(bucket_entry));
    bucket *buckets = malloc(256 * sizeof(bucket));
    if (!entries || !buckets) {
        free(entries);
        free(buckets);
        return;
    }
    init_buckets(buckets);
    maxlen = init_entries(strings, len, entries);
    head = &entries[0];
    /* Main loop */
    for (c = 0; c < maxlen-1; c++) {
        // printf("Bucketing on char %d from the right\n", c);
        bucket_strings(head, buckets, c);
        // print_buckets(buckets);
        head = merge_buckets(buckets);
        // print_list(head);
        init_buckets(buckets);
    }
    /* Copy back into array */
    copy_list(head, strings);
    free(buckets);
    free(entries);
}

void substring(char *s, char *sub, int p, int l) {
   int c = 0;
   
   while (c < l) {
      sub[c] = s[p+c-1];
      c++;
   }
   sub[c] = '\0';
}

void boss(int *C, int *last, char *W, int *Wi, char **kmers, int line, int k){
    int i, j;
    int pos = 0;
    char last_node[k-1];
    char last_Wi[k-2];
    int Wi_freq[255];

    strcpy(last_node, "");
    strcpy(last_Wi, "");
    memset(Wi_freq, 0,  sizeof(int)*255);

    for(i = 0; i < line; i++){
        if (strcmp(kmers[i], kmers[i+1]) != 0){
            int same_wi = 0;
            char label[k-1];
            substring(kmers[i], label, k-2, k-1);
            char label_Wi[k-2];
            substring(kmers[i], label_Wi, k-2, k-2);
            
            if(strcmp(last_Wi, label_Wi) != 0){
                memset(Wi_freq, 0, sizeof(int)*255);   
            } else {
                same_wi = 1;
            }   
            if(kmers[i][0] == '$'){
                int same_kmer = 1;
                for(j = 1; j < k; j++){
                    if(kmers[i][j] != kmers[i+1][j]){
                        same_kmer = 0;
                        break;
                    }
                }
                if(!same_kmer){
                    if(Wi_freq[kmers[i][0]] > 0){
                        Wi[pos] = 0;
                    } else {
                        Wi[pos] = 1;
                        Wi_freq[kmers[i][0]]++;
                    }
                    if(strcmp(label, last_node) == 0){
                        last[pos-1] = 0;
                        last[pos] = 1;
                    } else {
                        last[pos] = 1;
                    }
                    W[pos] = kmers[i][0];
                    C[kmers[i][0]]++;
                    pos++;
                    substring(kmers[i], last_node, k-2, k-2);
                    substring(kmers[i], last_node, k-2, k-1);
                }
            } else {
                if(Wi_freq[kmers[i][0]] > 0){
                    Wi[pos] = 0;
                } else {
                    Wi[pos] = 1;
                    Wi_freq[kmers[i][0]]++;
                }
                if(strcmp(label, last_node) == 0){
                    last[pos-1] = 0;
                    last[pos] = 1;
                } else {
                    last[pos] = 1;
                }
                W[pos] = kmers[i][0];
                C[kmers[i][0]]++;
                pos++;
                substring(kmers[i], last_Wi, k-2, k-2);
                substring(kmers[i], last_node, k-2, k-1);
            }
        }
    }

    C[1] = C['$'];
    C[2] = C['A'] + C[1];
    C[3] = C['C'] + C[2];
    C[4] = C['G'] + C[3];
    C[5] = C['N'] + C[4];

    C[0] = 0;

    char alphabet[6] = {'$', 'A', 'C', 'G', 'N', 'T'};
    printf("C array:\n");
    for(j = 0; j < 6; j++){
        printf("%c %d\n", alphabet[j], C[j]);
    }
    printf("\n");

    printf("\nBOSS:\n");
    printf("last\tW\t\tWi\n");
    for(i = 0; i < pos; i++){
        printf("%d\t\t%c\t\t%d\n", last[i], W[i], Wi[i]);
    }

}

int compare(const void *a, const void *b) {
    const char *c1 = *(char**)a;
    const char *c2 = *(char**)b;
    return c1[0] > c2[0];
}

int main(int argc, char *argv[]) {
    int i, j;

    if (argc < 3) {
        printf("Missing argument!\nusage: ./bgc [file.fastq] [k]\n");
        exit(-1);
    }

    FILE *fastqFile = fopen(argv[1], "r");
    FILE *reads = fopen("reads.txt", "w");
    int k = atoi(argv[2]);

    char *buff = (char*)malloc(105*sizeof(char));
    int qtd = 0;
    while(1) {
        fscanf(fastqFile, " %[^\n]", buff);
        fscanf(fastqFile, " %[^\n]", buff);
        // for(i = 0; i < k; i++){
        //     fprintf(reads, "$");
        // }
        fprintf(reads, "%s$\n", buff);
        fscanf(fastqFile, " %[^\n]", buff);
        fscanf(fastqFile, " %[^\n]", buff);
        qtd++;
        if(feof(fastqFile)) { 
            break ;
        }     
    }
    fclose(fastqFile);
    fclose(reads);

    reads = fopen("reads.txt", "r");

    fseek(reads, 0, SEEK_END);
    int n = ftell(reads);
    rewind(reads);

    k++;


    char **kmers = (char**)malloc(n*sizeof(char*));
    for(i = 0; i < n; i++)
        kmers[i] = (char*)malloc((k+1)*sizeof(char));

    int line = 0;
    while(1) {
        fscanf(reads, "%s", buff);

        if(feof(reads)) { 
            break ;
        } 

        int len = strlen(buff);
        char kmer[k+1];
        for(i = 0; i < len-k+1; i++){
            for(j = 0; j < k; j++){
                kmer[k-j-1] = buff[i+j];
            }
            kmer[k] = '\0';
            strcpy(kmers[line], kmer);
            kmers[line][k] = '\0';
            line++;
        }
        
    }

    // for(i = 0; i < line; i++){
    //     printf("%s\n", kmers[i]);
    // }

    // printf("\nordenei pela aresta\n");
    
    qsort(kmers, line, sizeof(*kmers), compare);

    // for(i = 0; i < line; i++){
    //     printf("%s\n", kmers[i]);
    // }

    radix_sort_string(kmers, line);

    // printf("\nordenei com radix\n");

    // for(i = 0; i < line; i++){
    //     printf("%s\n", kmers[i]);
    // }
    // printf("\n");

    char W[line+1];
    int last[line+1];
    int Wi[line+1];
    int C[255];
    memset(C, 0, sizeof(int)*255);

    boss(C, last, W, Wi, kmers, line, k);
    
}
