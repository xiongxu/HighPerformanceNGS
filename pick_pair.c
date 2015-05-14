// gcc -g -O3 -Wall pick_pair.c -o pick_pair  -I/share/software/software/zlib_1.2.8_install/include -L/share/software/software/zlib_1.2.8_install/lib -lz
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <zlib.h>
#include <time.h>
#include <sys/time.h>

#define myprintf_Fastq(line,out)	\
	do{								\
		if ((line)){				\
			gzprintf(out,"%s\n%s\n+\n%s",(line)->name,(line)->seq,(line)->quality);	\
			free((line)->name);		\
			free((line)->seq);		\
			free((line)->quality);	\
			free((line));			\
		}							\
	}while (0)

typedef struct _fastq_
{
	char *name;
	char *seq;
	char *quality;
}Fastq;

struct globalArgs_t {
	char *read1;
	char *read2;
	const char *outfile;
} globalArgs;

static inline Fastq *readNextNode(gzFile fq,char *buf);
gzFile creat_gzinfile(const char *infile);
gzFile creat_gzoutfile(const char *outfile,const char *suffix);
int load_fastq_file(const char *fq_file1,const char *fq_file2,const char *outfile,long long begin);
void display_usage(char * argv[]);
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
		*(buf+strlen(buf)-1)=0;	//remove the last \n in the buf
//		*(buf+(strchr(buf,' ')-buf))='\0';
		line->name  = strdup(buf);
		
		buf=gzgets(fq,buf,1024*sizeof(char));
		*(buf+strlen(buf)-1)=0;	//remove the last \n in the buf
		line->seq = strdup(buf);

		buf=gzgets(fq,buf,1024*sizeof(char));

		buf=gzgets(fq,buf,1024*sizeof(char));
		line->quality = strdup(buf);
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

int load_fastq_file(const char *fq_file1,const char *fq_file2,const char *outfile,long long begin)
{
	unsigned long total_reads_count[2];
	memset(total_reads_count,0,sizeof(unsigned long ));
	char *buf = (char *)malloc(1024 * sizeof(char));
	gzFile fastq1 = creat_gzinfile(fq_file1);
	gzFile fastq2 = creat_gzinfile(fq_file2);
	gzFile out1=creat_gzoutfile(outfile,"_1_PE.fq.gz");
	gzFile out2=creat_gzoutfile(outfile,"_1_SE.fq.gz");
	gzFile out3=creat_gzoutfile(outfile,"_2_PE.fq.gz");
	gzFile out4=creat_gzoutfile(outfile,"_2_SE.fq.gz");
	
	Fastq *line1=NULL,*line2=NULL;
	while (1) {
		line1=readNextNode(fastq1,buf);
		line2=readNextNode(fastq2,buf);
		while (line1 && strncmp(line1->name,line2->name,(strchr(line1->name,' ')-line1->name))<0){
			myprintf_Fastq(line1,out2);
			line1=readNextNode(fastq1,buf);
		}
		while (line2 && strncmp(line1->name,line2->name,(strchr(line1->name,' ')-line1->name))>0){
			myprintf_Fastq(line2,out4);
			line2=readNextNode(fastq2,buf);
		}
		if (line1==NULL && line2==NULL) break;
		myprintf_Fastq(line1,out1);
		myprintf_Fastq(line2,out3);
	}
	fprintf(stderr,"Finished load file at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	gzclose(fastq1);
	gzclose(fastq2);
	gzclose(out1);
	gzclose(out2);
	gzclose(out3);
	gzclose(out4);
	free(buf);
	return 0;
}

void display_usage(char * argv[]){
	char *buffer=(char* )malloc(10240*sizeof(char));
	const char* usage=
"\nCopyright (c) 2015\n" \
"Contact: XiongXu <xuxiong19880610@163.com> <xiongxu@me.com> \n" \
"Discription:\n  This program is used for picking out the paired fastq format reads.\n" \
"Usage: %s [-1 READ1] [-2 READ2] [-o OUTFILE] [-h] \n" \
"Example1:\n  %s -1 /home/litc/work/mate_pair/data/Maize/M017_L7_I030_R1.clean.fastq -2 /home/litc/work/mate_pair/data/Maize/M017_L7_I030_R2.clean.fastq -o pair_out \n" \
"\n" \
"   [-1 READ1]  = fastq formated file1.                                [required]\n" \
"   [-2 READ2]  = fastq formated file2.                                [option]\n" \
"   [-o OUTPUT] = OUTPUT file.                                         [required]\n" \
"   [-h]        = This helpful help screen.                            [option]\n" \
"\n";
	sprintf(buffer,usage,argv[0],argv[0]);
	fprintf(stderr,"%s",buffer);
	free(buffer);
	exit(1);
}

int main(int argc, char *argv[]){
	int opt = 0;
	globalArgs.read1=NULL;
	globalArgs.read2=NULL;
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
	load_fastq_file(globalArgs.read1,globalArgs.read2,globalArgs.outfile,begin);
	fprintf(stderr,"Finished  at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	return 0;
}
