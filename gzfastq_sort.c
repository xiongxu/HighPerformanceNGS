//gcc -g -O3 -Wall gzfastq_sort.c -I/share/software/software/zlib_1.2.8_install/include -L/share/software/software/zlib_1.2.8_install/lib -o gzfastq_sort -lz
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <getopt.h>
#include <err.h>
#include <time.h>
#include <sys/time.h>
#include "IO_stream.h"

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
	unsigned long reads_num;
} globalArgs;

void display_usage(char * argv[]);
unsigned long str2unsigned_long(char *str);
static inline int comp_seq(const void *a, const void *b);
static inline int comp_name(const void *a, const void *b);
void load_file(gzFile fp,const char *outfile,unsigned long *Line_num,unsigned long *max_read_num);
static inline Fastq *readNextNode(gzFile fq,char *buf);
void executeCMD(const char *cmd, char *result);
static long long usec(void);

long long usec(void){
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return (((long long)tv.tv_sec)*1000000)+tv.tv_usec;
}

void display_usage(char * argv[]){
	char *buffer=(char* )malloc(10240*sizeof(char));
	const char* usage=
"\nCopyright (c)  2015\n" \
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

unsigned long str2unsigned_long(char *str){
	unsigned long temp = 0;
	if (*str == '-' ) {
		fprintf(stderr,"reads count must be a positive integer!\n");
		exit(1);
	}
	while(*str != 0) {//*str!='\0'
		if ((*str < '0') || (*str > '9')) break;
		temp = temp * 10 + (*str - '0');
		str++;
	}
	return temp;
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

void load_file(gzFile fp,const char *outfile,unsigned long *Line_num,unsigned long *max_read_num) {
	Fastq *chunk=(Fastq *)calloc(*max_read_num,sizeof(Fastq));
	char *buf = (char *) malloc(1024*sizeof(char));
	Fastq *line=NULL;
	unsigned int i=0;
	long long begin=usec();
	while (1) {
		line=readNextNode(fp,buf);
		if (line==NULL) break;
//		fprintf(stderr,"sizeof(line): %d\tsizeof(Fastq): %d\tsizeof(Fastq *): %d\n",sizeof(line),sizeof(Fastq),sizeof(Fastq *));
		memmove(chunk+*Line_num,line,sizeof(Fastq));
		free(line);
		(*Line_num)++;
	}
	fprintf(stderr,"done read file at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	
	FILE *out=NULL;
	if (globalArgs.by_seq==0 && globalArgs.by_name==1) {
		qsort(chunk,*Line_num,sizeof(Fastq),comp_name);
		out=fcreat_outfile(outfile,"_sort_by_name.fq");
	}
	else if (globalArgs.by_seq==1 && globalArgs.by_name==0) {
		qsort(chunk,*Line_num,sizeof(Fastq),comp_seq);
		out=fcreat_outfile(outfile,"_sort_by_seq.fq");
	}

	fprintf(stderr,"done qsort file at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);

	for (i=0;i<*Line_num ;i++ ) {
		fprintf(out,"%s\n%s\n+\n%s\n",(chunk+i)->name,(chunk+i)->seq,(chunk+i)->quality);
		free_Fastq(chunk+i);
	}
	free(buf);
	free(chunk);
	fclose(out);
	fprintf(stderr,"done write file at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
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

void executeCMD(const char *cmd, char *result){
	char buf_ps[1024];
	char ps[1024]={0};
	FILE *ptr;
	strcpy(ps, cmd);
	if((ptr=popen(ps, "r"))!=NULL){
		while(fgets(buf_ps, 1024, ptr)!=NULL){
			strcat(result, buf_ps);
			if(strlen(result)>1024) break;
		}
		pclose(ptr);
		ptr = NULL;
	}
	else{
		printf("popen %s error\n", ps);
	}
}

unsigned long count_read(gzFile fq){
	char *buf = (char *)malloc(1024 * sizeof(char));
	unsigned long i=0;
	while (gzgets(fq,buf,1024*sizeof(char))!=NULL){
		gzgets(fq,buf,1024*sizeof(char));
		gzgets(fq,buf,1024*sizeof(char));
		gzgets(fq,buf,1024*sizeof(char));
		i++;
	}
	gzrewind(fq);
	free(buf);
	fprintf(stderr,"total_reads_num: %ld\n",i);
	return i;
}

int main(int argc, char *argv[]) {
	int opt = 0;
	globalArgs.infile="-";
	globalArgs.outfile="-";
	globalArgs.by_name=0;
	globalArgs.by_seq=0;
	globalArgs.reads_num=0;
	const char *optString = "i:o:r:nsh?";
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
			case 'r':
				globalArgs.reads_num = str2unsigned_long(optarg);
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
	unsigned long max_line_num;
	if (globalArgs.reads_num==0){
		max_line_num = count_read(fp);
		fprintf(stderr,"max_reads_num: %ld\n",max_line_num);
	}
	else{
		max_line_num = globalArgs.reads_num;
	}

	fprintf(stderr,"name: %d\tseq: %d\n",globalArgs.by_name,globalArgs.by_seq);
	unsigned long reads_num=0;
	load_file(fp,globalArgs.outfile,&reads_num,&max_line_num);
	gzclose(fp);
	return 0;
}
