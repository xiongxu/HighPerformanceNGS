//  gcc -g -O3 -Wall gzfastq_uniq_sort.c -o gzfastq_uniq_sort -I/share/software/software/zlib_1.2.8_install/include -I. -L/share/software/software/zlib_1.2.8_install/lib -L. -l hashtbl -lz
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <zlib.h>
#include <time.h>
#include <sys/time.h>
#include "hashtbl.h"

#define free_Name_Quality(line)	\
	do{							\
		free(line.name);		\
		free(line.quality);	\
	}while (0)

#define free_Fastq(line)	\
	do{							\
		free((line)->seq);							\
		free_Name_Quality((line)->Name_quality);	\
		free((line));								\
	}while (0)

typedef struct _Name_Quality_
{
	char *name;
	char *quality;
}Name_Quality;

typedef struct _fastq_
{
	Name_Quality Name_quality;
	char *seq;
}Fastq;

typedef struct _new_fastq_
{
	unsigned long count;
	Name_Quality uniq_fastq[2];
}New_Fastq;

struct globalArgs_t {
	const char *read1;
	const char *read2;
	short ends;
	const char *outfile;
} globalArgs;

static inline Fastq *readNextNode(gzFile fq,char *buf);
gzFile creat_gzinfile(const char *infile);
gzFile creat_gzoutfile(const char *outfile,const char *suffix);
HASHTBL *load_fastq_file(const char *fq_file1,const char *fq_file2,unsigned long *total_reads_count,long long begin,short *strLen);
unsigned long count_read(gzFile fq,char *buf);
void display_usage(char * argv[]);
static inline int compare_hashed_data_count(const void *a, const void *b);
void qsort_output_hash(HASHTBL *hashtbl,const char *outfile,short strLen);
void output_hash(HASHTBL *hashtbl,const char *outfile,short strLen);
static inline long long usec(void);

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
		*(buf+strlen(buf)-1)='\0';	//remove the last \n in the buf
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

gzFile creat_gzinfile(const char *infile) {
	gzFile fiz=gzopen(infile,"rb");
	if(fiz == NULL) {
		fprintf(stderr,"open file %s failed\n",infile);
		exit(1);
	}
	return fiz;
}

gzFile creat_gzoutfile(const char *outfile,const char *suffix) {
	char *out=(char *)calloc(strlen(outfile)+128,sizeof(char));
	sprintf(out,"%s%s",outfile,suffix);
	gzFile fo = gzopen(out,"wb");
	if (fo == NULL){
		fprintf(stderr,"open file %s failed\n",out);
		exit(1);
	}
	return fo;
}

HASHTBL *load_fastq_file(const char *fq_file1,const char *fq_file2,unsigned long *total_reads_count,long long begin,short *strLen)
{
	char *buf = (char *)malloc(1024 * sizeof(char));
	gzFile fastq1=creat_gzinfile(fq_file1);

	HSIZE elecnt=(HSIZE)count_read(fastq1,buf);
//	fprintf(stderr,"Finished %s count_read at %.3f s\n",fq_file1,(double)(usec()-begin)/CLOCKS_PER_SEC);
	HSIZE tblsiz=1.34*elecnt;
	gzFile fastq2 = NULL;
	if (globalArgs.ends>1){
		fastq2=creat_gzinfile(fq_file2);
	}
	HASHTBL *tbl = hashtbl_create(tblsiz, NULL);
	
	Fastq *line1=NULL,*line2=NULL;
	char *pair_seq=(char *)calloc(1024 , sizeof(char));
	while (1){
		line1=readNextNode(fastq1,buf);
		if (line1==NULL) break;
		if (! *strLen) *strLen=(short)strlen(line1->seq);
//		if (strlen(line1->seq)>75) {
//			strncpy(pair_seq,line1->seq,50);
//		}else{
			strncpy(pair_seq,line1->seq,strlen(line1->seq));
//		}
		if (globalArgs.ends>1) {
			line2=readNextNode(fastq2,buf);
			if (line2==NULL || strncmp(line1->Name_quality.name,line2->Name_quality.name,strchr(line1->Name_quality.name,' ')-line1->Name_quality.name)){
				fprintf(stderr,"error at %ld: %s\n",*total_reads_count,line1->Name_quality.name);
				break;
			}
//			if (strlen(line2->seq)>75) {
//				strncpy(pair_seq+50,line2->seq,50);
//			}else{
				strncpy(pair_seq+strlen(line1->seq),line2->seq,strlen(line2->seq));
//			}
		}
		void *data=hashtbl_get(tbl,pair_seq);
		if (data==NULL){
			New_Fastq *new_fastq=(New_Fastq *)malloc(sizeof(New_Fastq));
			new_fastq->uniq_fastq[0]=line1->Name_quality;
			if (globalArgs.ends>1) new_fastq->uniq_fastq[1]=line2->Name_quality;
			new_fastq->count=1;
			hashtbl_insert(tbl, pair_seq, (void *)new_fastq);
		}
		else{
			((*(New_Fastq *)data).count)++;
			free_Fastq(line1);
			if (globalArgs.ends>1) free_Fastq(line2);
		}
		(*total_reads_count)++;
		if(*total_reads_count % (elecnt/10)==0){
			fprintf(stderr,"loaded %lu at %.3f s\n",*total_reads_count,(double)(usec()-begin)/CLOCKS_PER_SEC);
		}
		memset(pair_seq,0,1024 * sizeof(char));
	}
	gzclose(fastq1);
	free(pair_seq);
	if (globalArgs.ends>1) gzclose(fastq2);
	fprintf(stderr,"unique reads number = %d\n",(int)tbl->count);
	fprintf(stderr,"Finished load hash at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	return tbl;
}

unsigned long count_read(gzFile fq,char *buf){
	unsigned long i=0;
	while (gzgets(fq,buf,1024*sizeof(char))!=NULL){
		gzgets(fq,buf,1024*sizeof(char));
		gzgets(fq,buf,1024*sizeof(char));
		gzgets(fq,buf,1024*sizeof(char));
		i++;
	}
	gzrewind(fq);
	fprintf(stderr,"total_reads_num: %ld\n",i);
	return i;
}

void display_usage(char * argv[]){
	char *buffer=(char* )malloc(10240*sizeof(char));
	const char* usage=
"\nCopyright (c) Berrygenomics 2012-2013\n" \
"Contact: XiongXu <xuxiong@berrygenomics.com> <xiongxu@me.com> \n" \
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

int compare_hashed_data_count(const void *a, const void *b) {
	return ((New_Fastq *)(*(ENTRY* const *)b)->data)->count-((New_Fastq *)(*(ENTRY* const *)a)->data)->count;
}

void qsort_output_hash(HASHTBL *hashtbl,const char *outfile,short strLen){
	gzFile out1=creat_gzoutfile(outfile,"_1_uniq.fq.gz");
	gzFile out2=NULL;
	if (globalArgs.ends>1) out2=creat_gzoutfile(outfile,"_2_uniq.fq.gz");
	ENTRY **nodes=dump_hash_table(hashtbl);
//	fprintf(stderr,"Start qsort\n");
	qsort(nodes, hashtbl->count, sizeof(ENTRY *), compare_hashed_data_count);
	unsigned long i;
	char Sequence[strLen+1];
	for (i=0;i<hashtbl->count ;i++ ) {
		New_Fastq *temp_data=(New_Fastq *)(nodes[i]->data);
		memcpy(Sequence,nodes[i]->key,strLen*sizeof(char));
		Sequence[strLen]=0;
		if (globalArgs.ends>1) {
			gzprintf(out2,"%s\t%ld\n%s\n+\n%s\n",temp_data->uniq_fastq[1].name,temp_data->count,
				nodes[i]->key+strLen,temp_data->uniq_fastq[1].quality);
			free_Name_Quality(temp_data->uniq_fastq[1]);
		}
		gzprintf(out1,"%s\t%ld\n%s\n+\n%s\n",temp_data->uniq_fastq[0].name,temp_data->count,
		Sequence,temp_data->uniq_fastq[0].quality);
		free_Name_Quality(temp_data->uniq_fastq[0]);
		free(temp_data);
	}
	free(nodes);
	gzclose(out1);
	if (globalArgs.ends>1) gzclose(out2);
}

void output_hash(HASHTBL *hashtbl,const char *outfile,short strLen){
	gzFile out1=creat_gzoutfile(outfile,"_1_uniq.fq.gz");
	gzFile out2=NULL;
	if (globalArgs.ends>1) out2=creat_gzoutfile(outfile,"_2_uniq.fq.gz");
	unsigned long i;
	ENTRY *entry;
	char Sequence[strLen+1];
	for (i = 0; i < hashtbl->size; ++i) {
		entry = hashtbl->nodes[i];
		while (entry) {
			New_Fastq *temp_data=(New_Fastq *)(entry->data);
			memcpy(Sequence,entry->key,strLen*sizeof(char));
			Sequence[strLen]=0;
			if (globalArgs.ends>1) {
				gzprintf(out2,"%s\t%ld\n%s\n+\n%s\n",temp_data->uniq_fastq[1].name,temp_data->count,
					entry->key+strLen,temp_data->uniq_fastq[1].quality);
				free_Name_Quality(temp_data->uniq_fastq[1]);
			}
			gzprintf(out1,"%s\t%ld\n%s\n+\n%s\n",temp_data->uniq_fastq[0].name,temp_data->count,
				Sequence,temp_data->uniq_fastq[0].quality);
			free_Name_Quality(temp_data->uniq_fastq[0]);
			free(temp_data);
			entry = entry->next;
		}
	}
	gzclose(out1);
	if (globalArgs.ends>1) gzclose(out2);
}

void hash_data_delete(void *DATA) {
	ENTRY *node=(ENTRY *)DATA;
//	free(node->data);
	free(node);
}

int main(int argc, char *argv[]){
	int opt = 0;
	globalArgs.read1=NULL;
	globalArgs.read2=NULL;
	globalArgs.ends=0;
	globalArgs.outfile="out";
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
				globalArgs.outfile = optarg;
				globalArgs.ends++;
				break;
			case '2':
				globalArgs.read2 = optarg;
				globalArgs.ends++;
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

	unsigned long total_reads_num=0;
	fprintf(stderr,"%s",globalArgs.read1);
	if (globalArgs.read2){
		fprintf(stderr,"\t%s\n",globalArgs.read2);
	}else{
		fprintf(stderr,"\n");
	}
	short SeqLen=0;
	HASHTBL *hash = load_fastq_file(globalArgs.read1,globalArgs.read2,&total_reads_num,begin,&SeqLen);
	fprintf(stderr,"hash size: %ld\ntotal reads = %ld\n",(unsigned long)hash->size,total_reads_num);
	fprintf(stderr,"unique reads percentage: %.3f%%\n",(float)hash->count/total_reads_num*100);
	
//	output_hash(hash,globalArgs.outfile,SeqLen);
	qsort_output_hash(hash,globalArgs.outfile,SeqLen);
//	fprintf(stderr,"Finished qsort_output_hash at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	hashtbl_destroy(hash,hash_data_delete);
	fprintf(stderr,"Finished  at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	return 0;
}
