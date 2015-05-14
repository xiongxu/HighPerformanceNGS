//gcc -g -O3 -Wall gzfastq_sample.0.2.c -o gzfastq_sample_0_2 -lz -I/share/work1/staff/xuxiong/software/backup/fastq-tools-0.6/src -L/share/work1/staff/xuxiong/software/backup/fastq-tools-0.6/src -lrng
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <fcntl.h>
#include <getopt.h>
#include <err.h>
#include <time.h>
#include <sys/time.h>
#include <inttypes.h>
#include <unistd.h>
#include <libgen.h>
#include "khash.h"
#include "common.h"
#include "rng.h"

#define STR(x) #x

#define freeNode(line)	\
do{								\
	if(line){					\
		free((line)->name);		\
		free((line)->seq);		\
		free((line)->quality);	\
		free((line));			\
	}							\
}while (0)

#define printNode(Line,stream,i)	\
do{																												\
	if (Line){																									\
		if (globalArgs.fastq){gzprintf(stream,"%s_%lu\n%s\n+\n%s",(Line)->name,i,(Line)->seq,(Line)->quality);}	\
		if (globalArgs.fasta){gzprintf(stream,">%s_%lu\n%s\n",Line->name,i,Line->seq);}							\
		freeNode(Line);																							\
	}																											\
}while (0)

#define MAKE_GZFILE(TYPE)			\
gzFile creat_gzoutfile_##TYPE (const char *prefix,TYPE n,const char *suffix,void(*f)(char *,const char *,TYPE,const char *)){	\
	char *out=(char *)calloc(strlen(prefix)+128,sizeof(char));			\
	f(out,prefix,n,suffix);												\
	gzFile fo = gzopen(out,STR(wb));									\
	if (fo == NULL){													\
		fprintf(stderr,STR(open file %s failed\n),out);					\
		exit(1);														\
	}																	\
	return fo;															\
}

#define MAKE_FILE(TYPE)			\
FILE *creat_outfile_##TYPE (const char *prefix,TYPE n,const char *suffix,void(*f)(char *,const char *,TYPE,const char *)){	\
	char *out=(char *)calloc(strlen(prefix)+128,sizeof(char));			\
	f(out,prefix,n,suffix);												\
	FILE *fo = fopen(out,STR(wb));										\
	if (fo == NULL){													\
		fprintf(stderr,STR(open file %s failed\n),out);					\
		exit(1);														\
	}																	\
	return fo;															\
}

typedef struct _fastq_
{
	char *name;
	char *seq;
	char *quality;
}Fastq;

struct globalArgs_t {
	char *read1;
	char *read2;
	short ends;
	const char *outfile;
	unsigned long reads_n;
	short fasta;
	short fastq;
} globalArgs;

static uint32_t g_subsam_seed = 0;
static double g_subsam_frac = -1.;

void display_usage(char * argv[]);
static inline void sprintUL(char *outfile,const char *prefix,unsigned long n,const char *suffix);
static inline void sprintDouble(char *outfile,const char *prefix,double n,const char *suffix);
static inline gzFile creat_gzinfile(const char *infile);
static inline int filter_reads(const char *read_ID);
void shuffle(rng_t* rng, unsigned long* xs, unsigned long n);
unsigned long* index_without_replacement(rng_t* rng, unsigned long n);
static inline int cmpul( const void* a, const void* b );
void get_number_from_file(char *infile1,char *infile2,unsigned long pick_out,long long begin) ;
void proportion_file(char *infile1,char *infile2);
static inline Fastq *readNextNode(gzFile fq,char *buf) ;
void executeCMD(const char *cmd, char *result);
unsigned long get_reads_number(const char *infile);
unsigned long count_read(gzFile fq,char *buf);
static inline long long usec(void);

MAKE_GZFILE(double)
MAKE_GZFILE(khint64_t)
MAKE_FILE(double)
MAKE_FILE(khint64_t)

long long usec(void){
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return (((long long)tv.tv_sec)*1000000)+tv.tv_usec;
}

void display_usage(char * argv[]){
	char *buffer=(char* )malloc(10240*sizeof(char));
	const char* usage=
"\nCopyright (c) 2015\n" \
"Contact: XiongXu <xuxiong19880610@163.com> <xiongxu@me.com> \n" \
"Discription:\n  This program is used for get subsample from large text based gzip fastq files .\n" \
"Usage: %s {-1 fastq1} [-2 fastq2] [-o OUTFILE] [-s FLOAT] [-n UL] [-h] \n" \
"Example1:\n  %s  -1 /share/work1/staff/xuxiong/ChIP-seq/ChIP-Input_L4_I044.R1.clean.fastq.gz -2 /share/work1/staff/xuxiong/ChIP-seq/ChIP-Input_L4_I044.R2.clean.fastq.gz -s 1.01 \n" \
"\n" \
"   [-1 fastq1] = fastq1.                                                           [required]\n" \
"   [-2 fastq2] = fastq2.                                                           [option]\n" \
"   [-s FLOAT]  = fraction of templates to subsample; integer part as seed [-1].    [option]\n" \
"   [-n UL]     = number of picked reads,not work with -s.                          [option]\n" \
"   [-f ]       = output fasta format.                                              [option]\n" \
"   [-q ]       = output fastq format[default].                                     [option]\n" \
"   [-h]        = This helpful help screen.                                         [option]\n" \
"\n";
	sprintf(buffer,usage,argv[0],argv[0]);
	fprintf(stderr,"%s",buffer);
	free(buffer);
	exit(1);
}

void sprintUL(char *outfile,const char *prefix,unsigned long n,const char *suffix){
	sprintf(outfile,"%s.%lu%s",prefix,n,suffix);
}

void sprintDouble(char *outfile,const char *prefix,double n,const char *suffix){
	sprintf(outfile,"%s.%f%s",prefix,n,suffix);
}

gzFile creat_gzinfile(const char *infile) {
	gzFile fiz=gzopen(infile,"rb");
	if(fiz == NULL) {
		fprintf(stderr,"open file %s failed\n",infile);
		exit(1);
	}
	return fiz;
}

int filter_reads(const char *read_ID){
	uint32_t k = __ac_X31_hash_string(read_ID) + g_subsam_seed;
	return ((double)(k&0xffffff)/0x1000000 >= g_subsam_frac) ? 1:0;
}

/* randomly shuffle an array of unsigned longs */

void shuffle(rng_t* rng, unsigned long* xs, unsigned long n){
	unsigned long i, j, k;
	for (i = n - 1; i > 0; --i) {
		j = fastq_rng_uniform_int(rng, i + 1);
		k = xs[j]; xs[j] = xs[i]; xs[i] = k; // swap
	}
}

unsigned long* index_without_replacement(rng_t* rng, unsigned long n){
	unsigned long* xs = malloc_or_die(n * sizeof(unsigned long));
	size_t i;
	for (i = 0; i < n; ++i) xs[i] = i;
	shuffle(rng, xs, n);
	return xs;
}

int cmpul( const void* a, const void* b ) {
	return *(unsigned long*) a - *(unsigned long*) b ;
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

unsigned long get_reads_number(const char *infile){
	char *command=(char *)calloc(128,sizeof(char));
	fprintf(stderr,"%s\n",infile+(strlen(infile)-3));
	if (!strncmp(infile+(strlen(infile)-3),".gz",3)){
		sprintf(command,"zegrep -c '^[ATCGN]+$' %s",infile);
	}
	else{
		sprintf(command,"egrep -c '^[ATCGN]+$' %s",infile);
	}
	fprintf(stderr,"%s\n",command);
	char *char_max_line_num=(char *)calloc(128,sizeof(char));
	executeCMD(command,char_max_line_num);
	unsigned long max_line_num= strtoul(char_max_line_num,NULL,10);
	fprintf(stderr,"total_reads_num: %ld\n",max_line_num);
	free(char_max_line_num);
	free(command);
	return max_line_num;
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

void get_number_from_file(char *infile1,char *infile2,unsigned long pick_out,long long begin) {
	char *buf=(char *)calloc(1024,sizeof(char));
	gzFile fiz1=creat_gzinfile(infile1);
	gzFile fout1=creat_gzoutfile_khint64_t(basename(infile1),pick_out,".gz",sprintUL);
	unsigned long n=count_read(fiz1,buf);
//	unsigned long n=get_reads_number(globalArgs.read1);
	fprintf(stderr,"Finished count_read at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	if (pick_out > n) {
		fprintf(stderr,"pick_count > read_count (%lu > %lu)\n",pick_out,n);
		exit(0);
//		pick_out = reads_count;
	}
	gzFile fiz2=NULL,fout2=NULL;
	if (infile2!=NULL) {
		fiz2=creat_gzinfile(infile2);
		fout2=creat_gzoutfile_khint64_t(basename(infile2),pick_out,".gz",sprintUL);
	}
	Fastq *Line1=NULL,*Line2=NULL;
	unsigned long i=0,j=0,rng_seed = 4357;
	
	rng_t* rng = fastq_rng_alloc();
	fastq_rng_seed(rng, rng_seed);
	unsigned long *xs=index_without_replacement(rng, n);
	qsort(xs, pick_out, sizeof(unsigned long), cmpul);
	fprintf(stderr,"Start_read at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	while (1){
		Line1=readNextNode(fiz1,buf);
		if (j >= pick_out || Line1 == NULL) break;
		if (globalArgs.ends>1) Line2=readNextNode(fiz2,buf);
		if (j < pick_out && xs[j] == i) {
			printNode(Line1,fout1,i+1);
			printNode(Line2,fout2,i+1);
			++j;
		}
		else{
			freeNode(Line1);
			freeNode(Line2);
		}
		i++;
	}
	fprintf(stderr,"End_read at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	free(buf);
	gzclose(fout1);
	gzclose(fiz1);
	if (infile2!=NULL) {
		gzclose(fiz2);
		gzclose(fout2);
	}
	fastq_rng_free(rng);
	free(xs);
	fprintf(stderr,"total reads: %lu\npick out: %lu (%lu/%lu=%.6f)\n",n,pick_out,pick_out,n,(double)pick_out/n);
}

void proportion_file(char *infile1,char *infile2) {
	gzFile fiz1=creat_gzinfile(infile1);
	gzFile fout1=creat_gzoutfile_double(basename(infile1),g_subsam_frac,".gz",sprintDouble);
	gzFile fiz2=NULL,fout2=NULL;
	if (infile2!=NULL) {
		fiz2=creat_gzinfile(infile2);
		fout2=creat_gzoutfile_double(basename(infile2),g_subsam_frac,".gz",sprintDouble);
	}
	char *buf=(char *)calloc(1024,sizeof(char));
	Fastq *Line1=NULL,*Line2=NULL;
	unsigned long n=0,pick_out=0;
	while (1){
		Line1=readNextNode(fiz1,buf);
		if (Line1 == NULL) break;
		n++;
		if (globalArgs.ends>1) Line2=readNextNode(fiz2,buf);
		if (filter_reads(Line1->name)){
			freeNode(Line1);
			freeNode(Line2);
		}else{
			pick_out++;
			printNode(Line1,fout1,n);
			printNode(Line2,fout2,n);
		}
	}
	free(buf);
	gzclose(fout1);
	gzclose(fiz1);
	if (infile2!=NULL) {
		gzclose(fiz2);
		gzclose(fout2);
	}
	fprintf(stderr,"total reads: %lu\npick out: %lu (%lu/%lu=%.6f)\n",n,pick_out,pick_out,n,(double)pick_out/n);
}

Fastq *readNextNode(gzFile fq,char *buf) {
	Fastq *line = (Fastq *)malloc(sizeof(Fastq));
	buf=gzgets(fq,buf,1024*sizeof(char));
	if (!gzeof(fq)) {
		*(buf+strlen(buf)-1)='\0';
		line->name  = strdup(buf);
		
		buf=gzgets(fq,buf,1024*sizeof(char));
		*(buf+strlen(buf)-1)='\0';	//remove the last \n in the buf
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

int main(int argc, char *argv[]) {
	int opt = 0;
	globalArgs.read1=NULL;
	globalArgs.read2=NULL;
	globalArgs.outfile="out";
	globalArgs.ends=0;
	globalArgs.reads_n=0;
	globalArgs.fastq=1;
	globalArgs.fasta=0;
	const char *optString = "1:2:o:s:n:qfh?";
	char *q;
	if (argc<2) display_usage(argv);
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
		switch( opt ) {
			case '1':
				globalArgs.read1 = optarg;
				globalArgs.ends++;
				break;
			case '2':
				globalArgs.read2 = optarg;
				globalArgs.ends++;
				break;
			case 'o':
				globalArgs.outfile = optarg;
				break;
			case 's':
				if ((g_subsam_seed = strtol(optarg, &q, 10)) != 0) {
					srand(g_subsam_seed);
					g_subsam_seed = rand();
				}
				g_subsam_frac = strtod(q, &q);
				break;
			case 'n':
				globalArgs.reads_n = strtoul(optarg, NULL, 10);
				break;
			case 'f':
				globalArgs.fasta++;
				globalArgs.fastq=0;
				break;
			case 'q':
				globalArgs.fastq++;
				globalArgs.fasta=0;
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

	if (g_subsam_frac>0){
		proportion_file(globalArgs.read1,globalArgs.read2);
	}
	if (globalArgs.reads_n){
		get_number_from_file(globalArgs.read1,globalArgs.read2,globalArgs.reads_n,begin);
	}
	fprintf(stderr,"Finished at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	return 0;
}
