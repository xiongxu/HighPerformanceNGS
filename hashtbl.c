// gcc -g -O3 -Wall  -c hashtbl.c && ar -rcs libhashtbl.a hashtbl.o && rm hashtbl.o
#include <string.h>
#include <stdio.h>
#include "hashtbl.h"


/*
static HSIZE bkdr_hash(const char *key)
{
	HSIZE seed = 131; // 31 131 1313 13131 131313 etc.
	HSIZE hash = 0;

	while (*key)
		hash = hash * seed + (unsigned char)(*key++);

	return hash;
}

static inline HSIZE __ac_X31_hash_string(const char *s)
{
	HSIZE h = (UBYTE)*s;
	if (h) {
		for (++s ; *s; ++s) 
			h = (h << 5) - h + (UBYTE)*s;
	}
	return h;
}
*/
static HSIZE dictGenHashFunction(const char *buf) {
	size_t len=strlen(buf);
	HSIZE hash = 5381;
	while (len--)
		hash = ((hash << 5) + hash) + (unsigned char)(*buf++); /* hash * 33 + c */
	return hash;
}

HASHTBL *hashtbl_create(HSIZE size, HSIZE (*hashfunc)(const char *))
{
	HASHTBL *hashtbl;
	if (!(hashtbl=malloc(sizeof(HASHTBL))))
		return NULL;
	if (!(hashtbl->nodes = calloc(size, sizeof(ENTRY*)))) {
		free(hashtbl);
		return NULL;
	}
	hashtbl->size = size;
	hashtbl->count = 0;
	if (hashfunc){
		hashtbl->hashfunc = hashfunc;
	}else{
		hashtbl->hashfunc = dictGenHashFunction;
	}
	return hashtbl;
}
/*
void hashtbl_destroy(HASHTBL *hashtbl)
{
	HSIZE n;
	ENTRY *node, *oldnode;

	for (n = 0; n < hashtbl->size; ++n) {
		node = hashtbl->nodes[n];
		while (node) {
			free(node->key);
			oldnode = node;
			node = node->next;
			free(oldnode);
		}
	}
	free(hashtbl->nodes);
	free(hashtbl);
}
*/
void hashtbl_destroy(HASHTBL *hashtbl,void (*DELETE)(void *))
{
	HSIZE n;
	ENTRY *node, *oldnode;
	if (DELETE) {
		for (n = 0; n < hashtbl->size; ++n) {
			node = hashtbl->nodes[n];
			while (node!=NULL) {
				free(node->key);
				oldnode = node;
				node = node->next;
				hashtbl->count--;
				DELETE(oldnode);
			}
		}
	}else{
		for (n = 0; n < hashtbl->size; ++n) {
			node = hashtbl->nodes[n];
			while (node!=NULL) {
				free(node->key);
				oldnode = node;
				node = node->next;
				hashtbl->count--;
				free(oldnode);
			}
		}
	}
	free(hashtbl->nodes);
	free(hashtbl);
}

int hashtbl_insert(HASHTBL *hashtbl, const char *key, void *data)
{
	/* resize table if the threshold is exceeded
	 * default threshold is:
	 * <table size> * <load factor 0.75> */
	if (hashtbl->count >= hashtbl->size * 0.75) {
		hashtbl_resize(hashtbl, hashtbl->size * 2 + 1);
	}
	ENTRY *node;
	HSIZE hash = hashtbl->hashfunc(key) % hashtbl->size;
	node = hashtbl->nodes[hash];
	/* check if the key is already in the hashtbl */
	while (node) {
		if (!strcmp(node->key, key)) {
			node->data = data;
//			fprintf(stderr,"exist node->key: %s\n",node->key);
			return 0;
		}
		node = node->next;
	}
	/* create new entry */
	if (!(node = malloc(sizeof(ENTRY)))){
		fprintf(stderr,"malloc node error\n");
		return -1;
	}
	if (!(node->key = strdup(key))) {
		fprintf(stderr,"malloc node->key error\n");
		free(node);
		return -1;
	}
	node->data = data;
	node->next = hashtbl->nodes[hash];
	hashtbl->nodes[hash] = node;
	hashtbl->count++;
	return 0;
}

int hashtbl_remove(HASHTBL *hashtbl, const char *key)
{
	ENTRY *node, *prevnode = NULL;
	HSIZE hash = hashtbl->hashfunc(key) % hashtbl->size;

	node = hashtbl->nodes[hash];
	while (node) {
		if (!strcmp(node->key, key)) {
			free(node->key);
			if (prevnode){
				prevnode->next = node->next;
			}else{
				hashtbl->nodes[hash] = node->next;
			}
			free(node);
			hashtbl->count--;
			return 0;
		}
		prevnode = node;
		node = node->next;
	}

	return -1;
}

void *hashtbl_get(HASHTBL *hashtbl, const char *key)
{
	if (!hashtbl || !hashtbl->size) return NULL;
	ENTRY *node;
	HSIZE hash = hashtbl->hashfunc(key) % hashtbl->size;
	node = hashtbl->nodes[hash];
	while (node) {
		if (!strcmp(node->key, key)) {
			return node->data;
		}
		node = node->next;
	}
	return NULL;
}

int hashtbl_resize(HASHTBL *hashtbl, HSIZE size)
{
	HASHTBL newtbl;
	HSIZE n;
	ENTRY *node;
	fprintf(stderr,"resize at: %d\t%d\n",(int)hashtbl->size,(int)size);
	newtbl.size = size;
	newtbl.count = 0;
	newtbl.hashfunc = hashtbl->hashfunc;

	if (!(newtbl.nodes = calloc(size, sizeof(ENTRY*)))){
		fprintf(stderr,"malloc node error\n");
		return -1;
	}

	for (n = 0; n < hashtbl->size; ++n) {
		for (node = hashtbl->nodes[n]; node; node = node->next) {
			hashtbl_insert(&newtbl, node->key, node->data);
			hashtbl_remove(hashtbl, node->key);
		}
	}

	free(hashtbl->nodes);
	hashtbl->size = newtbl.size;
	hashtbl->count = newtbl.count;
	hashtbl->nodes = newtbl.nodes;

	return 0;
}

void hashtbl_dump(HASHTBL *hashtbl, FILE *stream)
{
	unsigned int i;
	ENTRY *entry;
	for (i = 0; i < hashtbl->size; ++i) {
		entry = hashtbl->nodes[i];
		if (entry)
			fprintf(stream,"\n[%u]", i);
		while (entry) {
			fprintf(stream," -> (\"%s\":%p)", entry->key, entry->data);
			entry = entry->next;
		}
	}
	printf("\n");
}

ENTRY** dump_hash_table(HASHTBL *hashtbl)
{
	ENTRY** D = (ENTRY**)malloc(hashtbl->count * sizeof(ENTRY *));
	unsigned int i,j;
	ENTRY *entry;
	for (i = 0,j=0; i < hashtbl->size; ++i) {
		entry = hashtbl->nodes[i];
		while (entry) {
			D[j++]=entry;
			entry = entry->next;
		}
	}
	return D;
}

UBYTE *int2char(int pos,UBYTE *value){
	int i;
	for (i=3;i>=0 ;i-- )
		value[i]= (pos>>(i*7) & 0x7F) +1;
	value[4]=0;
	return value;
}

int char2int(UBYTE *string,int *pos){
	int i;
	*pos=0;
	for (i=3;i>=0 ;i-- ){
		*pos<<=7;
		*pos+=(string[i]-1);
	}
	return *pos;
}


int cmp_char2int(const void *a, const void *b){
    return atoi(*(char * const*)a)-atoi(*(char * const*)b);
}

int cmp_twoBitInt(const void *a, const void *b){
    int c=0,d=0;
    return char2int(*(UBYTE * const*)a,&c)-char2int(*(UBYTE * const*)b,&d);
}

int compare_hashed_data_key(const void *a, const void *b) {
    return strcmp((*(ENTRY* const *)a)->key,(*(ENTRY* const *)b)->key);
}

char **union_hashed_keys(HASHTBL *hashtbl_a,HASHTBL *hashtbl_b){
	char **keys=(char **)calloc((hashtbl_a->count+hashtbl_b->count),sizeof(char *));
	unsigned int i,j=0;
	ENTRY *entry;
	for (i = 0; i < hashtbl_a->size; ++i) {
		entry = hashtbl_a->nodes[i];
		while (entry) {
			keys[j++]=entry->key;
			entry = entry->next;
		}
	}
	for (i = 0; i < hashtbl_b->size; ++i) {
		entry = hashtbl_b->nodes[i];
		while (entry) {
			keys[j++]=entry->key;
			entry = entry->next;
		}
	}
//	fprintf(stderr,"keys_count: %d\t%d\n",j,hashtbl_a->count+hashtbl_b->count);
//	qsort(keys,j,sizeof(char *),cmp_char2int);
	qsort(keys,j,sizeof(char *),cmp_twoBitInt);
	return keys;
}
