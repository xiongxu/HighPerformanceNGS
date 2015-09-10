//gcc -shared -std=gnu99 -fPIC -g -O2 -DNDEBUG convertUniqueBin2Pos.c -o convertUniqueBin2Pos.so  -I/share/software/software/R-3.0_install/R-3.0.1/include -I/usr/local/include -I/home/xuxiong/work/HighPerformanceNGS/zlib-1.2.8/include -I/home/xuxiong/work/HighPerformanceNGS/hiredis -L/usr/local/lib64 -L/share/software/software/R-3.0_install/R-3.0.1/lib -L/home/xuxiong/work/HighPerformanceNGS/zlib-1.2.8/lib -L/home/xuxiong/work/HighPerformanceNGS/hiredis -lR -lz -lhiredis
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <getopt.h>
#include <zlib.h>
#include <time.h>
#include <sys/time.h>
#include "sds.h"
#include "IO_stream.h"
#include "split.h"
#include "dict.h"
#include "list.h"

#define free_Name_Quality(line)	\
do{							\
	free((line).name);		\
	free((line).quality);	\
}while (0)

#define free_Fastq(line)	\
do{												\
	free((line)->seq);							\
	free_Name_Quality((line)->Name_quality);	\
	free((line));								\
}while (0)

#define initdictTypeSE(mydictType)	\
do{													\
	(mydictType)->hashFunction=myhashFunction;		\
	(mydictType)->keyDup=mykeyDup;					\
	(mydictType)->valDup=NULL;						\
	(mydictType)->keyCompare=mykeyCompare;			\
	(mydictType)->keyDestructor=mykeyDestructor;	\
	(mydictType)->valDestructor=myvalDestructorSE;	\
}while (0)

#define initdictTypePE(mydictType)	\
do{													\
	(mydictType)->hashFunction=myhashFunction;		\
	(mydictType)->keyDup=mykeyDup;					\
	(mydictType)->valDup=NULL;						\
	(mydictType)->keyCompare=mykeyCompare;			\
	(mydictType)->keyDestructor=mykeyDestructor;	\
	(mydictType)->valDestructor=myvalDestructorPE;	\
}while (0)

#define SUMQuality(sumQ,string,strLen)	\
do{											\
	uint32_t index;							\
	for (index=0;index<strLen ;index++ ){	\
		sumQ+=(uint8_t) *(string+index);	\
	}										\
}while (0)

#define AllocDataSE(line1,Count,SumQuality)	\
Data_t *data=(Data_t *)malloc(sizeof(Data_t));	\
do{												\
	data->uniq_fastq=list_create(freeNameQuality, comp_name );		\
	list_add_data(data->uniq_fastq,&line1->Name_quality);					\
	data->count=Count;							\
	data->sumQuality=SumQuality;				\
}while (0)

#define printNameQualityList(chunk,stream)	\
do{													\
	Node *n = chunk->head;							\
	while (n) {										\
		Node *next_n = n->next;						\
		Name_Quality *NQ=(Name_Quality *)n->data;	\
		fprintf(stream,"%s\n",NQ->quality);			\
		n = next_n;									\
	}												\
}while(0)

#define printSortedDict(ht,suffix,Outfile)   \   
do{                                                                 \
    if (ht){                                                        \
        dictEntry** D=dump_dict(ht);                                \
        FILE *stream=fcreat_outfile(Outfile,suffix);                \
        unsigned long i=0;                                          \
        for (i=0;i<ht->used;i++){                                   \
            Data_t *value=(Data_t *)(D[i]->val);					\
			fprintf(stream,"%s\t%u\n%s\n+\n",((Name_Quality *)value->uniq_fastq->head->data)->name,value->count,(char *)D[i]->key);	\
			printNameQualityList(value->uniq_fastq,stream);						\
        }                                                           \
        free(D);                                                    \
        fclose(stream);                                             \
    }                                                               \
}while(0)

typedef struct _Name_Quality_
{
	char *name;
	char *quality;
	uint16_t seq_len;
}Name_Quality;

typedef struct _fastq_
{
	Name_Quality Name_quality;
	char *seq;
}Fastq;

typedef struct _new_fastq_
{
	uint32_t sumQuality;
	uint32_t count;
	List *uniq_fastq;
}Data_t;

struct globalArgs_t {
	const char *read1;
	int sortbyseq;
	int sortbyCount;
	const char *outfile;
} globalArgs;

static inline Fastq *readNextNode(gzFile fq,char *buf);
dict *load_fastq_SE(const char *fq_file1,unsigned long *total_reads_count,dictType *mydictType,long long begin);

void display_usage(char * argv[]);
static inline long long usec(void);

static inline void *mykeyDup(void *privdata, const void *key);
static inline int mykeyCompare(void *privdata, const void *key1, const void *key2);
static inline void mykeyDestructor(void *privdata, void *key);
static inline void myvalDestructorSE(void *privdata, void *obj);

static inline unsigned int myhashFunction(const void *key);
static inline int compare_hashed_key(const void *a, const void *b);
static inline int compare_hashed_data_count(const void *a, const void *b);
static inline void freeNameQuality(Name_Quality *line);
static inline int comp_name(const void *a, const void *b);

unsigned int myhashFunction(const void *key){
	return dictGenHashFunction((unsigned char *)key,sdslen((sds) key));
}

void *mykeyDup(void *privdata, const void *key){
	return (void *)sdsdup((sds) key);
}

int mykeyCompare(void *privdata, const void *key1, const void *key2){
	return sdscmp((sds) key1,(sds) key2);
}

void mykeyDestructor(void *privdata, void *key){
	sdsfree((sds) key);
}

void myvalDestructorSE(void *privdata, void *obj){
	Data_t *myobj=(Data_t *)obj;
	list_free(myobj->uniq_fastq);
	free(myobj);
}

void freeNameQuality(Name_Quality *line){
	free_Name_Quality(*line);
	free(line);
}

long long usec(void){
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return (((long long)tv.tv_sec)*1000000)+tv.tv_usec;
}

int comp_name(const void *a, const void *b) {
	Name_Quality *c=(Name_Quality *)a;
	Name_Quality *d=(Name_Quality *)b;
	if(strlen(c->name)!=strlen(d->name)) {
		return strlen(c->name)-strlen(d->name);
	}else{
		return strcmp(c->name, d->name);
	}
}

Fastq *readNextNode(gzFile fq,char *buf) {
	Fastq *line = (Fastq *)malloc(sizeof(Fastq));
	buf=gzgets(fq,buf,1024*sizeof(char));
	if (!gzeof(fq)) {
		*(buf+strlen(buf)-1)='\0';
		line->Name_quality.name  = strdup(buf);
		
		buf=gzgets(fq,buf,1024*sizeof(char));
		line->Name_quality.seq_len=strlen(buf)-1;
		*(buf+line->Name_quality.seq_len)='\0';	//remove the last \n in the buf
		line->seq = strdup(buf);

		buf=gzgets(fq,buf,1024*sizeof(char));

		buf=gzgets(fq,buf,1024*sizeof(char));
		*(buf+strlen(buf)-1)='\0';
		line->Name_quality.quality = strdup(buf);
		return line;
	} else {
		free(line);
		return NULL;
	}
}

dict *load_fastq_SE(const char *fq_file1,unsigned long *total_reads_count,dictType *mydictType,long long begin){
	char *buf = (char *)malloc(1024 * sizeof(char));
	gzFile fastq1=open_input_stream(fq_file1);
	dict *ht=dictCreate(mydictType,NULL);
	Fastq *line1=NULL;
	sds pair_seq=sdsempty();
	uint32_t sumQ;
	while (1){
		line1=readNextNode(fastq1,buf);
		if (line1==NULL) break;
		pair_seq=sdsnewlen(line1->seq,line1->Name_quality.seq_len);
		sumQ=0;
		SUMQuality(sumQ,line1->Name_quality.quality,line1->Name_quality.seq_len);
		dictEntry *entry=dictFind(ht,pair_seq);
		if (entry==NULL){
			AllocDataSE(line1,1,sumQ);
			dictAdd(ht, pair_seq,data);
			free(line1->seq);
		}else{
			Data_t *value=(Data_t *)entry->val;
			value->count++;
			if (sumQ > value->sumQuality){
				value->sumQuality=sumQ;
			}
			list_add_data(value->uniq_fastq, &line1->Name_quality);
			free(line1->seq);
		}
		sdsfree(pair_seq);
		(*total_reads_count)++;
	}
	gzclose(fastq1);
	fprintf(stderr,"unique reads number = %lu(%lu / %lu = %.3f%%)\n",ht->used,ht->used,*total_reads_count,100.0*ht->used/(*total_reads_count));
	fprintf(stderr,"hash size: %ld\n",ht->size);
	fprintf(stderr,"Finished load hash at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	return ht;
}

int compare_hashed_data_count(const void *a, const void *b) {
    return ((Data_t *)(*(dictEntry* const *)b)->val)->count-((Data_t *)(*(dictEntry* const *)a)->val)->count;
}

int compare_hashed_key(const void *a, const void *b) {
	return sdscmp((sds)(*(dictEntry* const *)a)->key,(sds)(*(dictEntry* const *)b)->key);
}

dictEntry** dump_dict(dict *ht){
    dictEntry** D = (dictEntry **)malloc(ht->used * sizeof(dictEntry *));
    dictIterator *iter=dictGetIterator(ht);
    dictEntry *entry=dictNext(iter);
    uint32_t j=0;
    for (; entry; entry=dictNext(iter)) {
        D[j++]=entry;
    }
    dictReleaseIterator(iter);
    qsort(D, ht->used, sizeof(dictEntry *),(globalArgs.sortbyseq & !globalArgs.sortbyCount ? compare_hashed_key : compare_hashed_data_count));
    return D;
}

void display_usage(char * argv[]){
	char *buffer=(char* )malloc(10240*sizeof(char));
	const char* usage=
"\nCopyright (c) 2015\n" \
"Contact: XiongXu <xuxiong19880610@163.com> <xiongxu@me.com> \n" \
"Discription:\n  If you haven't appointed the read2 file, this program just handle the read1 file. Otherwise, it will also handle the read2 fastq file, and treat each of the joined read1 and read2 sequence as the key, 52193025 with 100 nt length pair end reads will cost about 22G memery(with 44%% unique pair end reads).\n" \
"Usage: %s [-1 READ1] [-C sort by seq count] [-S sort by seq][-o OUTFILE] [-h] \n" \
"Example1:\n  %s -1 /home/litc/work/mate_pair/data/Maize/M017_L7_I030_R1.clean.fastq -S -o pair_out \n" \
"Example2:\n  %s -1 /home/litc/work/mate_pair/data/Maize/M017_L7_I030_R1.clean.fastq -C -o single_out \n" \
"\n" \
"   [-1 READ1]  = fastq formated file1.                                [required]\n" \
"   [-C ]       = sort by sequence count.                              [option]\n" \
"   [-S ]       = sort by sequence.                                    [option]\n" \
"   [-o OUTPUT] = OUTPUT file.                                         [required]\n" \
"   [-h]        = This helpful help screen.                            [option]\n" \
"\n";
	sprintf(buffer,usage,argv[0],argv[0],argv[0]);
	fprintf(stderr,"%s",buffer);
	free(buffer);
	exit(1);
}

int main(int argc, char *argv[]){
	int opt = 0;
	globalArgs.read1="-";
	globalArgs.sortbyseq=1;
	globalArgs.sortbyCount=0;
	globalArgs.outfile="-";
	const char *optString = "1:o:CSh?";
	if (argc<2) {
		display_usage(argv);
		exit(1);
	}
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
		switch( opt ) {
			case '1':
				globalArgs.read1 = optarg;
				break;
			case 'S':
				globalArgs.sortbyseq = 1;
				globalArgs.sortbyCount = 0;
				break;
			case 'C':
				globalArgs.sortbyseq = 0;
				globalArgs.sortbyCount = 1;
				break;
			case 'o':
				globalArgs.outfile = optarg;
				break;
			case '?':	/* fall-through is intentional */
			case 'h':
				display_usage(argv);
				break;
			default:
				fprintf(stderr,"error parameter!\n");
				break;
		}
		opt = getopt( argc, argv, optString );
	}
	long long begin;
	begin=usec();

	dictType *mydictType=(dictType *)malloc(sizeof(dictType));
	dict *hash=NULL;

	unsigned long total_reads_num=0;
	initdictTypeSE(mydictType);
	hash = load_fastq_SE(globalArgs.read1,&total_reads_num,mydictType,begin);
	printSortedDict(hash,"_sortKeyUniq.fq",globalArgs.outfile);
	dictRelease(hash);
	free(mydictType);
	fprintf(stderr,"Finished  at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	return 0;
}
