// gcc -g -O3 -Wall fastq_trim.c -o fastq_trim -I/share/software/software/zlib_1.2.8_install/include -L/share/software/software/zlib_1.2.8_install/lib -lz
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <sys/time.h>
#include <time.h>
#include "IO_stream.h"

#define free_Fastq(line)	\
	do{							\
		free((line)->name);		\
		free((line)->seq);		\
		free((line)->quality);	\
	}while (0)

struct globalArgs_t {
	const char *infile;	/* -i option */
	const char *outfile;	/* -o option */
	int start;
	int end;
	int verbosity;	/* -v option */
	int gzip;
} globalArgs;

typedef struct _fastq_
{
	char *name;
	char *seq;
	char *quality;
}Fastq;

void load_file(gzFile fp,FILE *out,unsigned long *Line_num,int S,int E);
static inline int readNextNode(gzFile fq,Fastq *line,char *buf,int S,int E) ;
void display_usage(char * argv[]);
static inline long long usec(void);

long long usec(void){
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return (((long long)tv.tv_sec)*1000000)+tv.tv_usec;
}

void display_usage(char * argv[]){
	char *buffer=(char* )malloc(10240*sizeof(char));
	const char* usage=
"\nCopyright (c) 2015\n" \
"Contact: XiongXu <xuxiong@19880610@163.cn> <xiongxu@me.com> \n" \
"Discription:\n  This program is used for cutting the reads to get the specified cycles.\n" \
"Usage: %s [-i Infile] [-o OUTFILE] [-s start] [-e end] [-h] \n" \
"Example1:\n  %s -i /share/work1/staff/xuxiong/test/13C37198_L7_I012.R1.clean.fastq.gz -s 0 -e 100 -o out \n" \
"Example2:\n  zcat /share/work1/staff/xuxiong/test/13C37198_L7_I012.R1.clean.fastq.gz |%s -e 100 -o sss \n" \
"Example3:\n  zcat /share/work1/staff/xuxiong/test/13C37198_L7_I012.R1.clean.fastq.gz |%s -e 50 |less -S  \n" \
"\n" \
"   [-i Infile]    = Infile.default is stdin                        [required]\n" \
"   [-o OUTPUT]    = OUTPUT file.default is stdout.                 [option]\n" \
"   [-s Start]     = 0 based start position,default is 0            [option]\n" \
"   [-e End]       = 1 based end position,default is 400            [option]\n" \
"   [-h] This helpful help screen.                                  [option]\n" \
"\n";
	sprintf(buffer,usage,argv[0],argv[0],argv[0],argv[0]);
	fprintf(stderr,"%s",buffer);
	free(buffer);
	exit(1);
}

int readNextNode(gzFile fq,Fastq *line,char *buf,int S,int E) {
	int i=1;
	buf=gzgets(fq,buf,1024*sizeof(char));
	if (!gzeof(fq)) {
		*(buf+strlen(buf)-1)='\0';
		line->name  =  strdup(buf);
		
		buf=gzgets(fq,buf,1024*sizeof(char));
		*(buf+strlen(buf)-1)='\0';	//remove the last \n in the buf
		line->seq =calloc(E-S+1,sizeof(char));
		strncpy(line->seq,buf+S,E-S);

		buf=gzgets(fq,buf,1024*sizeof(char));
		
		buf=gzgets(fq,buf,1024*sizeof(char));
		*(buf+strlen(buf)-1)='\0';	//remove the last \n in the buf
		line->quality =calloc(E-S+1,sizeof(char));
		strncpy(line->quality,buf+S,E-S);
	} else {
		return 0;
	}
	return i;
}

void load_file(gzFile fp,FILE *out,unsigned long *Line_num,int S,int E) {
	long long begin,end;
	begin=usec();
	Fastq *line = (Fastq *)malloc(sizeof(Fastq));
	char *buf = (char *) malloc(1024*sizeof(char));
	while (1) {
		memset(buf,0,1024*sizeof(char));
		int status=readNextNode(fp,line,buf,S,E);
		if (!status) break;
		(*Line_num)++;
		fprintf(out,"%s\n%s\n+\n%s\n",line->name,line->seq,line->quality);
		free_Fastq(line);
	}
	end=usec();
	fprintf(stderr,"Total_reads: %lu\nFinished in %.3f s\n",*Line_num,(double)(end-begin)/CLOCKS_PER_SEC);
	free(buf);
	free(line);
}

int main(int argc, char *argv[])
{
	int opt = 0;
	globalArgs.infile="-";
	globalArgs.outfile="-";
	globalArgs.start=0;
	globalArgs.end=400;
	globalArgs.verbosity=0;
	globalArgs.gzip=0;
	const char *optString = "i:o:s:e:vzh?";
	if (argc<2){
		display_usage(argv);
		exit(1);
	}
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
		switch( opt ) {
			case 'i':
				globalArgs.infile = optarg;
				break;
			case 'o':
				globalArgs.outfile = optarg;
				break;
			case 'v':
				globalArgs.verbosity++;
				break;
			case 'z':
				globalArgs.gzip++;
				break;
			case 's':
				globalArgs.start = atoi(optarg);
				break;
			case 'e':
				globalArgs.end = atoi(optarg);
				break;
			case '?':	/* fall-through is intentional */
			case 'h':
				display_usage(argv);
				break;
			default:
				fprintf(stderr,"error parameter!\n");
				/* You won't actually get here. */
				break;
		}
		opt = getopt( argc, argv, optString );
	}
	gzFile in=open_input_stream(globalArgs.infile);
	FILE *out=fcreat_outfile(globalArgs.outfile,".trim.fastq");
	unsigned long reads_num=0;
	load_file(in,out,&reads_num,globalArgs.start,globalArgs.end);
	gzclose(in);
	fclose(out);
	return 0;
}
