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

#define AllocDataPE(line1,line2,Count,SumQuality)	\
Data_t *data=(Data_t *)malloc(sizeof(Data_t));							\
do{																		\
	data->uniq_fastq=(Name_Quality *)malloc(2*sizeof(Name_Quality));	\
	data->uniq_fastq[0]=line1->Name_quality;							\
	data->uniq_fastq[1]=line2->Name_quality;							\
	data->count=Count;													\
	data->sumQuality=SumQuality;										\
}while (0)

#define AllocDataSE(line1,Count,SumQuality)	\
Data_t *data=(Data_t *)malloc(sizeof(Data_t));	\
do{												\
	data->uniq_fastq=&line1->Name_quality;		\
	data->count=Count;							\
	data->sumQuality=SumQuality;				\
}while (0)

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
	Name_Quality *uniq_fastq;
}Data_t;

struct globalArgs_t {
	const char *read1;
	const char *read2;
	const char *outfile;
} globalArgs;

static inline Fastq *readNextNode(gzFile fq,char *buf);
void output_hashSE(dict *hashtbl,const char *outfile);
void output_hashPE(dict *hashtbl,const char *outfile);
dict *load_fastq_SE(const char *fq_file1,unsigned long *total_reads_count,dictType *mydictType,long long begin);
dict *load_fastq_PE(const char *fq_file1,const char *fq_file2,unsigned long *total_reads_count,dictType *mydictType,long long begin);

void display_usage(char * argv[]);
static inline long long usec(void);

static inline void *mykeyDup(void *privdata, const void *key);
static inline int mykeyCompare(void *privdata, const void *key1, const void *key2);
static inline void mykeyDestructor(void *privdata, void *key);
static inline void myvalDestructorSE(void *privdata, void *obj);
static inline void myvalDestructorPE(void *privdata, void *obj);
static inline unsigned int myhashFunction(const void *key);

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

void myvalDestructorPE(void *privdata, void *obj){
	Data_t *myobj=(Data_t *)obj;
	free_Name_Quality(myobj->uniq_fastq[0]);
	free_Name_Quality(myobj->uniq_fastq[1]);
	free(myobj->uniq_fastq);
	free(myobj);
}

void myvalDestructorSE(void *privdata, void *obj){
	Data_t *myobj=(Data_t *)obj;
	free_Name_Quality(*(myobj->uniq_fastq));
	free(myobj);
}

long long usec(void){
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return (((long long)tv.tv_sec)*1000000)+tv.tv_usec;
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

dict *load_fastq_PE(const char *fq_file1,const char *fq_file2,unsigned long *total_reads_count,dictType *mydictType,long long begin){
	char *buf = (char *)malloc(1024 * sizeof(char));
	gzFile fastq1=open_input_stream(fq_file1);
	gzFile fastq2=open_input_stream(fq_file2);
	
	dict *ht=dictCreate(mydictType,NULL);
	Fastq *line1=NULL,*line2=NULL;
	sds pair_seq=sdsempty();
	uint32_t sumQ;
	while (1){
		line1=readNextNode(fastq1,buf);
		if (line1==NULL) break;
		line2=readNextNode(fastq2,buf);
		if (line2==NULL || strncmp(line1->Name_quality.name,
			line2->Name_quality.name,strchr(line1->Name_quality.name,' ')-line1->Name_quality.name) ){
			fprintf(stderr,"error at %ld: %s\nunmatched read name\n",*total_reads_count,line1->Name_quality.name);
			break;
		}
		pair_seq=sdsnewlen(line1->seq,line1->Name_quality.seq_len);
		pair_seq=sdscatlen(pair_seq,line2->seq,line2->Name_quality.seq_len);
		sumQ=0;
		SUMQuality(sumQ,line1->Name_quality.quality,line1->Name_quality.seq_len);
		SUMQuality(sumQ,line2->Name_quality.quality,line2->Name_quality.seq_len);
		dictEntry *entry=dictFind(ht,pair_seq);
		if (entry==NULL){
			AllocDataPE(line1,line2,1,sumQ);
			dictAdd(ht, pair_seq,data);
			free(line1->seq);
			free(line2->seq);
		}else{
			Data_t *value=(Data_t *)entry->val;
			value->count++;
			if (sumQ > value->sumQuality){
				AllocDataPE(line1,line2,value->count,sumQ);
				dictReplace(ht, entry->key, data);
			}else{
				free_Fastq(line1);
				free_Fastq(line2);
			}
		}
		sdsfree(pair_seq);
		(*total_reads_count)++;
	}
	gzclose(fastq1);
	gzclose(fastq2);
	fprintf(stderr,"unique reads number = %lu(%lu / %lu = %.3f%%)\n",ht->used,ht->used,*total_reads_count,100.0*ht->used/(*total_reads_count));
	fprintf(stderr,"hash size: %ld\n",ht->size);
	fprintf(stderr,"Finished load hash at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	return ht;
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
				AllocDataSE(line1,value->count,sumQ);
				dictReplace(ht, entry->key, data);
			}else{
				free_Fastq(line1);
			}
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

void display_usage(char * argv[]){
	char *buffer=(char* )malloc(10240*sizeof(char));
	const char* usage=
"\nCopyright (c) 2015\n" \
"Contact: XiongXu <xuxiong19880610@163.com> <xiongxu@me.com> \n" \
"Discription:\n  If you haven't appointed the read2 file, this program just handle the read1 file. Otherwise, it will also handle the read2 fastq file, and treat each of the joined read1 and read2 sequence as the key, 52193025 with 100 nt length pair end reads will cost about 22G memery(with 44%% unique pair end reads).\n" \
"Usage: %s [-1 READ1] [-2 READ2] [-o OUTFILE] [-h] \n" \
"Example1:\n  %s -1 /home/litc/work/mate_pair/data/Maize/M017_L7_I030_R1.clean.fastq -2 /home/litc/work/mate_pair/data/Maize/M017_L7_I030_R2.clean.fastq -o pair_out \n" \
"Example2:\n  %s -1 /home/litc/work/mate_pair/data/Maize/M017_L7_I030_R1.clean.fastq -o single_out \n" \
"\n" \
"   [-1 READ1]  = fastq formated file1.                                [required]\n" \
"   [-2 READ2]  = fastq formated file2.                                [option]\n" \
"   [-o OUTPUT] = OUTPUT file.                                         [required]\n" \
"   [-h]        = This helpful help screen.                            [option]\n" \
"\n";
	sprintf(buffer,usage,argv[0],argv[0],argv[0]);
	fprintf(stderr,"%s",buffer);
	free(buffer);
	exit(1);
}

void output_hashPE(dict *hashtbl,const char *outfile){
//	gzFile out1=creat_outfile(outfile,"_1_uniq.fq.gz");
	FILE *out1=fcreat_outfile(outfile,"_1_uniq.fq");
//	gzFile out2=creat_outfile(outfile,"_2_uniq.fq.gz");
	FILE *out2=fcreat_outfile(outfile,"_2_uniq.fq");
	dictIterator *iter=dictGetIterator(hashtbl);
	dictEntry *entry=dictNext(iter);
	for (;entry;entry=dictNext(iter) ){
		Data_t *temp_data=(Data_t *)(entry->val);
		char **splitSeq=(char **)malloc(2*sizeof(char *));
		int string_index=split(splitSeq,temp_data->uniq_fastq[0].name,' ');
		sds key1=sdsnewlen((sds)entry->key,temp_data->uniq_fastq[0].seq_len);
//		gzprintf(out1,"%s XI:Z:%s\tXF:i:%u\n%s\n+\n%s\n",splitSeq[0],splitSeq[1],temp_data->count,key1,temp_data->uniq_fastq[0].quality);
//		fprintf(out1,"%s XI:Z:%s\tXF:i:%u\n%s\n+\n%s\n",splitSeq[0],splitSeq[1],temp_data->count,key1,temp_data->uniq_fastq[0].quality);
		fprintf(out1,"%s\t%u\n%s\n+\n%s\n",temp_data->uniq_fastq[0].name,temp_data->count,key1,temp_data->uniq_fastq[0].quality);
		sdsfree(key1);
		freeSplitString(splitSeq,string_index);

		splitSeq=(char **)malloc(2*sizeof(char *));
		string_index=split(splitSeq,temp_data->uniq_fastq[1].name,' ');
		sds key2=sdsnewlen((sds)entry->key+temp_data->uniq_fastq[0].seq_len,temp_data->uniq_fastq[1].seq_len);
//		gzprintf(out2,"%s XI:Z:%s\tXF:i:%u\n%s\n+\n%s\n",splitSeq[0],splitSeq[1],temp_data->count,key2,temp_data->uniq_fastq[1].quality);
//		fprintf(out2,"%s XI:Z:%s\tXF:i:%u\n%s\n+\n%s\n",splitSeq[0],splitSeq[1],temp_data->count,key2,temp_data->uniq_fastq[1].quality);
		fprintf(out2,"%s\t%u\n%s\n+\n%s\n",temp_data->uniq_fastq[1].name,temp_data->count,key2,temp_data->uniq_fastq[1].quality);
		sdsfree(key2);
		freeSplitString(splitSeq,string_index);
	}
	dictReleaseIterator(iter);
//	gzclose(out1);
//	gzclose(out2);
	fclose(out1);
	fclose(out2);
}

void output_hashSE(dict *hashtbl,const char *outfile){
//	gzFile out1=creat_outfile(outfile,"_uniq.fq.gz");
	FILE *out1=fcreat_outfile(outfile,"_uniq.fq");
	dictIterator *iter=dictGetIterator(hashtbl);
	dictEntry *entry=dictNext(iter);
	for (;entry;entry=dictNext(iter) ){
		Data_t *temp_data=(Data_t *)(entry->val);
		char **splitSeq=(char **)malloc(2*sizeof(char *));
		int string_index=split(splitSeq,temp_data->uniq_fastq->name,' ');
//		gzprintf(out1,"%s XI:Z:%s\tXF:i:%u\n%s\n+\n%s\n",splitSeq[0],splitSeq[1],temp_data->count,(char *)entry->key,temp_data->uniq_fastq->quality);
//		fprintf(out1,"%s XI:Z:%s\tXF:i:%u\n%s\n+\n%s\n",splitSeq[0],splitSeq[1],temp_data->count,(char *)entry->key,temp_data->uniq_fastq->quality);
		fprintf(out1,"%s\t%u\n%s\n+\n%s\n",temp_data->uniq_fastq->name,temp_data->count,(char *)entry->key,temp_data->uniq_fastq->quality);
		freeSplitString(splitSeq,string_index);
	}
	dictReleaseIterator(iter);
//	gzclose(out1);
	fclose(out1);
}

int main(int argc, char *argv[]){
	int opt = 0;
	globalArgs.read1="-";
	globalArgs.read2=NULL;
	globalArgs.outfile="-";
	const char *optString = "1:2:o:h?";
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
			case '2':
				globalArgs.read2 = optarg;
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
	if (globalArgs.read2){
		initdictTypePE(mydictType);
		hash = load_fastq_PE(globalArgs.read1,globalArgs.read2,&total_reads_num,mydictType,begin);
		output_hashPE(hash,globalArgs.outfile);
	}else{
		initdictTypeSE(mydictType);
		hash = load_fastq_SE(globalArgs.read1,&total_reads_num,mydictType,begin);
		output_hashSE(hash,globalArgs.outfile);
	}
	dictRelease(hash);
	free(mydictType);
	fprintf(stderr,"Finished  at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	return 0;
}
