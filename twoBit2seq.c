//gcc -g -O3 -Wall twoBit2seq.c -I/share/software/software/zlib_1.2.8_install/include -L/share/software/software/zlib_1.2.8_install/lib -o twoBit2seq -lz
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <assert.h>
#include <sys/time.h>
#include "IO_stream.h"
#include "twoBit.h"

typedef struct _fastq_ {
	char *name;
	char *seq;
	char *third;
	char *quality;
} Fastq;

struct globalArgs_t {
	char *infile;
	char *outfile;
	int compress;
} globalArgs;

void display_usage(char * argv[]);
void load_file(FILE *fp,FILE *fo,unsigned long *Line_num,int compress_level);
static long long usec(void);

int ntValNoN[256];
char valToNt[5];
char revNt[256];

long long usec(void){
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return (((long long)tv.tv_sec)*1000000)+tv.tv_usec;
}

void display_usage(char * argv[]){
	char *buffer=(char* )malloc(10240*sizeof(char));
	const char* usage=
"\nCopyright (c)  2012-2013\n" \
"Contact: XiongXu <xuxiong19880610@163.com> <xiongxu@me.com> \n" \
"Discription:\n  This program is used for decompressing and unpacking the gzip compressed 2bit file into ATCG sequence.\n" \
"Usage: %s [-i Infile] [-o OUTFILE] [-c compress_level][-h] \n" \
"Example1:\n  %s -i /share/work1/staff/xuxiong/test/13C37198_L7_I012.R1.clean.fastq.gz_sort_by_seq.2bit.gz -o 13C37198_L7_I012.R1.clean.fastq.gz_sort_by_seq.2bit.gz\n" \
"\n" \
"   [-i Infile] = Infile.default is stdin                     [required]\n" \
"   [-o OUTPUT] = OUTPUT file.default is stdout.              [option]\n" \
"   [-h] This helpful help screen.                            [option]\n" \
"\n";
	sprintf(buffer,usage,argv[0],argv[0]);
	fprintf(stderr,"%s",buffer);
	free(buffer);
	exit(1);
}

void load_file(FILE *fp,FILE *fo,unsigned long *Line_num,int compress_level) {
	char *buf = (char *)calloc(1024,sizeof(char));
	long long begin=usec();
	uint8_t seqlen=0,packedLen=0;
	char *unpacked_seq=NULL;
	fread(&seqlen,sizeof(uint8_t),1,fp);
	fread(&packedLen,sizeof(uint8_t),1,fp);
	while (1) {
		fread(buf,sizeof(UBYTE),packedLen,fp);
		if (feof(fp)) break;
		unpacked_seq=sds2seq((sds)buf,0,seqlen);
		fprintf(fo,"%s\n",unpacked_seq);
		free(unpacked_seq);
		(*Line_num)++;
	}
	fprintf(stderr,"done read file at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	free(buf);
}

int main(int argc, char *argv[]) {
	int opt = 0;
	globalArgs.infile="-";
	globalArgs.outfile="out";
	globalArgs.compress=4;
	const char *optString = "i:o:c:h?";
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
			case 'c':
				globalArgs.compress = atoi(optarg);
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
	initNtVal();
	FILE *fp=fopen_input_stream(globalArgs.infile);
	FILE *fo=fcreat_outfile(globalArgs.outfile,".decompress");
	unsigned long reads_num=0;
	load_file(fp,fo,&reads_num,globalArgs.compress);
	fclose(fp);
	fclose(fo);
	return 0;
}
