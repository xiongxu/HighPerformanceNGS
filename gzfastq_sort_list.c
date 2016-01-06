//gcc -g -O3 -Wall gzfastq_sort_list.c  -o gzfastq_sort_list -I/share/software/software/zlib_1.2.8_install/include -I. -L. -L/share/software/software/zlib_1.2.8_install/lib -lz -llist
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <getopt.h>
#include <err.h>
#include <time.h>
#include <sys/time.h>
#include "IO_stream.h"
#include "list.h"

typedef struct _fastq_ {
	char *name;
	char *seq;
	char *quality;
} Fastq;

#define free_Fastq(line)	\
	do{							\
		free((line)->name);		\
		free((line)->seq);		\
		free((line)->quality);	\
	}while (0)

struct globalArgs_t {
	const char *infile;
	const char *outfile;
	short by_name;
	short by_seq;
} globalArgs;

void display_usage(char * argv[]);
static inline void FastqFree(void *data);
static inline int comp_name(const void *a, const void *b);
static inline int comp_seq(const void *a, const void *b);
void load_file(gzFile fp,const char *outfile);
static inline Fastq *readNextNode(gzFile fq,char *buf);
static inline Fastq *dump_array(List *chunk,long long begin);
static long long usec(void);

long long usec(void){
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return (((long long)tv.tv_sec)*1000000)+tv.tv_usec;
}

void display_usage(char * argv[]){
	char *buffer=(char* )malloc(10240*sizeof(char));
	const char* usage=
"\nCopyright (c)  20155555\n" \
"Contact: XiongXu <xuxiong19880610@163.com> <xiongxu@me.com> \n" \
"Discription:\n  This program is used for sorting fastq files by name.\n" \
"Usage: %s [-i Infile] [-o OUTFILE] [-r reads_count_for allocating memory] [-s|-n] [-h] \n" \
"Example1:\n  zcat 13C37198_L7_I012.R1.clean.fastq.gz | %s -r 8000000 -o out -s\n" \
"\n" \
"   [-i Infile] = Infile.                                      [required]\n" \
"   [-o OUTPUT] = OUTPUT file.default is \'out\'.                 [option]\n" \
"   [-r reads_num] = reads_count_for allocating memory.         [option]\n" \
"   [-s ] sort by seqyence.                                     [option]\n" \
"   [-n ] sort by sequence name.                                [option]\n" \
"   [-h] This helpful help screen.                              [option]\n" \
"\n";
	sprintf(buffer,usage,argv[0],argv[0]);
	fprintf(stderr,"%s",buffer);
	free(buffer);
	exit(1);
}

void FastqFree(void *data){
	Fastq *line=(Fastq *)data;
	free_Fastq(line);
	free(line);
}

int comp_name(const void *a, const void *b) {
	Fastq *c=(Fastq *)a;
	Fastq *d=(Fastq *)b;
	if(strlen(c->name)!=strlen(d->name)) {
		return strlen(c->name)-strlen(d->name);
	}else{
		return strcmp(c->name, d->name);
	}
}

int comp_seq(const void *a, const void *b) {
	Fastq *c=(Fastq *)a;
	Fastq *d=(Fastq *)b;
	if(strlen(c->seq)!=strlen(d->seq)) {
		return strlen(c->seq)-strlen(d->seq);
	}else{
		return strcmp(c->seq, d->seq);
	}
}

Fastq *dump_array(List *chunk,long long begin){
	Fastq *arrays=(Fastq *)calloc(chunk->count,sizeof(Fastq));
	int i=chunk->count;
	Node *n = chunk->head;
	while (n) {
		Node *next_n = n->next;
		arrays[--i]=*(Fastq *)n->data;
		n = next_n;
	}
	fprintf(stderr,"done dump_array at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	qsort(arrays,chunk->count,sizeof(Fastq),(globalArgs.by_name && !globalArgs.by_seq) ? comp_name : comp_seq);
	return arrays;
}

void load_file(gzFile fp,const char *outfile) {
	List *chunk=list_create(FastqFree,(globalArgs.by_name && !globalArgs.by_seq) ? comp_name : comp_seq);
	char *buf = (char *) malloc(1024*sizeof(char));
	Fastq *line=NULL;
	long long begin=usec();
	while (1) {
		line=readNextNode(fp,buf);
		if (line==NULL) break;
		list_add_data(chunk,line);
	}
	fprintf(stderr,"done read file at %.3f s\nlist count: %d\n",(double)(usec()-begin)/CLOCKS_PER_SEC,chunk->count);
	
	FILE *out=(globalArgs.by_name && !globalArgs.by_seq) ? fcreat_outfile(outfile,"_sort_by_name.fq") : fcreat_outfile(outfile,"_sort_by_seq.fq");
	//list_sort(chunk);
	Fastq *arrays=dump_array(chunk,begin);
	fprintf(stderr,"done sort file at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);

//	Node *n = chunk->head;
//	while (n) {
//		Node *next_n = n->next;
//		line=(Fastq *)n->data;
//		fprintf(out,"%s\n%s\n+\n%s\n",line->name,line->seq,line->quality);
//		chunk->node_data_free_callback(n->data);
//		free(n);
//		n = next_n;
//	}
//	free(chunk);

	// while ( (line=(Fastq *)shift(chunk))!=NULL ){
	// 	fprintf(out,"%s\n%s\n+\n%s\n",line->name,line->seq,line->quality);
	// }
	int i;
	for (i=0;i<chunk->count ;i++ ) {
		fprintf(out,"%s\n%s\n+\n%s\n",(arrays+i)->name,(arrays+i)->seq,(arrays+i)->quality);
	}
	free(arrays);
	fprintf(stderr,"done write file at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	list_free(chunk);

	free(buf);
	fclose(out);
	fprintf(stderr,"done free list at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
}

Fastq *readNextNode(gzFile fq,char *buf) {
	Fastq *line = (Fastq *)malloc(sizeof(Fastq));
	buf=gzgets(fq,buf,1024*sizeof(char));
	if (!gzeof(fq)) {
		*(buf+strlen(buf)-1)=0;	//remove the last \n in the buf
//		*(buf+(strchr(buf,' ')-buf))=0;
		line->name  = strdup(buf);
		
		buf=gzgets(fq,buf,1024*sizeof(char));
		*(buf+strlen(buf)-1)=0;	//remove the last \n in the buf
		line->seq = strdup(buf);

		buf=gzgets(fq,buf,1024*sizeof(char));

		buf=gzgets(fq,buf,1024*sizeof(char));
		*(buf+strlen(buf)-1)=0;	//remove the last \n in the buf
		line->quality = strdup(buf);
		return line;
	} else {
		free(line);
		return NULL;
	}
}

int main(int argc, char *argv[]) {
	int opt = 0;
	globalArgs.infile="-";
	globalArgs.outfile="-";
	globalArgs.by_name=0;
	globalArgs.by_seq=0;
	const char *optString = "i:o:nsh?";
	if (argc<2) display_usage(argv);
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
		switch( opt ) {
			case 'i':
				globalArgs.infile = optarg;
				break;
			case 'o':
				globalArgs.outfile = optarg;
				break;
			case 'n':
				globalArgs.by_name = 1;
				globalArgs.by_seq = 0;
				break;
			case 's':
				globalArgs.by_name = 0;
				globalArgs.by_seq = 1;
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
	if (!globalArgs.by_name && !globalArgs.by_seq) globalArgs.by_seq=1;

	gzFile fp=open_input_stream(globalArgs.infile);
	fprintf(stderr,"name: %d\tseq: %d\n",globalArgs.by_name,globalArgs.by_seq);
	load_file(fp,globalArgs.outfile);
	gzclose(fp);
	return 0;
}
