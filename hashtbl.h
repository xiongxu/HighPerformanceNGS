#ifndef _HASHTBL_H
#define _HASHTBL_H

#include <stdio.h>
#include<stdlib.h>

//无论我们包含哪个标准头文件，<stddef.h>都会被自动包含进来
#define UBYTE unsigned char   /* Wants to be unsigned 8 bits. */
typedef size_t HSIZE;

typedef struct hashnode {
    char *key;
    void *data;
    struct hashnode *next;
} ENTRY;

typedef struct hashtbl {
    HSIZE size;
    HSIZE count;
    ENTRY **nodes;
    HSIZE (*hashfunc)(const char *);
} HASHTBL;


HASHTBL *hashtbl_create(HSIZE size, HSIZE (*hashfunc)(const char *));
void hashtbl_destroy(HASHTBL *hashtbl,void (*DELETE)(void *));
int hashtbl_insert(HASHTBL *hashtbl, const char *key, void *data);
int hashtbl_remove(HASHTBL *hashtbl, const char *key);
void *hashtbl_get(HASHTBL *hashtbl, const char *key);
int hashtbl_resize(HASHTBL *hashtbl, HSIZE size);
void hashtbl_dump(HASHTBL *hashtbl, FILE *stream);
ENTRY** dump_hash_table(HASHTBL *hashtbl);
char **union_hashed_keys(HASHTBL *hashtbl_a,HASHTBL *hashtbl_b);
int char2int(UBYTE *string,int *pos);
UBYTE *int2char(int pos,UBYTE *value);
int cmp_char2int(const void *a, const void *b);
int cmp_twoBitInt(const void *a, const void *b);
int compare_hashed_data_key(const void *a, const void *b);


#define DEF_SIZE 100000
#define HASHTBL_INIT(hashtbl) (hashtbl = hashtbl_create(DEF_SIZE, NULL))

#endif
