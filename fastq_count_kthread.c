//gcc -g -O3 -Wall klib-master/kthread.c fastq_count_kthread.c -o fastq_count_kthread -I/share/software/software/zlib_1.2.8_install/include -L/share/software/software/zlib_1.2.8_install/lib -lz -lpthread
#include "IO_stream.h"
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h>

#define CPU_NUM sysconf(_SC_NPROCESSORS_ONLN)

#define Alloc2DArray(Array,nrow,ncol)	\
do{																	\
	Array=(uint64_t **)calloc(nrow,sizeof(uint64_t *));				\
	uint32_t index=0;												\
	for (;index < nrow;index++ ){									\
		*(Array+index)=(uint64_t *)calloc(ncol,sizeof(uint64_t));	\
	}																\
}while (0)

#define Free2DArray(Array,nrow)	\
do{											\
	uint32_t index;							\
	for (index=0;index < nrow;index++ ){	\
		free(*(Array+index));				\
	}										\
	free(Array);							\
}while (0)

#define AssignQuality(Array,Q,Len)	\
do{													\
	uint32_t index;									\
	for(index=0;index<Len;index++){					\
		(*(*(Array+(uint8_t)*(Q+index))+index))++;	\
	}												\
}while (0)

#define statQ(Array,nrow,ncol,sum,Q1,sumQ1,Q2,sumQ2)	\
do{														\
	uint32_t ii=0,jj=0;									\
	for(ii=0;ii<nrow;ii++){								\
		for (jj=0;jj<ncol ;jj++ ){						\
			sum += *(*(Array+ii)+jj);					\
			if (ii >=Q1) sumQ1 += *(*(Array+ii)+jj);	\
			if (ii >=Q2) sumQ2 += *(*(Array+ii)+jj);	\
		}												\
	}													\
}while (0)

#define printSeqLen(Array,n,stream,minLen,maxLen,meanLen)	\
do{													\
	uint32_t index;									\
	fprintf(stream,"#Len:");						\
	for (index=minLen;index <= maxLen;index++ ){	\
		fprintf(stream,"\t%u",index);				\
	}												\
	fprintf(stream,"\n#Freq:");						\
	for (index=minLen;index <= maxLen;index++ ){	\
		fprintf(stream,"\t%lu",*(Array+index));		\
	}												\
	fprintf(stream,"\n");							\
}while (0)

#define statSeqLen(Array,n,minLen,maxLen,meanLen,sumFreq)	\
do{													\
	uint32_t index;									\
	for (index=0;index < n;index++ ){				\
		if (*(Array+index)){						\
			sumFreq+=*(Array+index);				\
			(meanLen)+=1.0*(*(Array+index))*index;	\
			if (!minLen) minLen=index;				\
			if (maxLen<index) maxLen=index;			\
		}											\
	}												\
}while (0)

struct globalArgs_t {
	char **infiles;
	const char *outfile;
	int numInfiles;
	int thread;
	int header;
	int LengthDetail;
} globalArgs;

typedef struct _thread_arg_
{
	char **infiles;
	FILE *out;
	uint32_t *read_count;
	double *mean_len;
}thread_arg_t;

void load_fq(void *data,long i,int tid);
static inline void display_usage(char * argv[]);
void count_read(gzFile fq,uint32_t *i,const char *infile,FILE *out,double *mean_length);
static inline long long usec(void);
void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);

long long usec(void){
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return (((long long)tv.tv_sec)*1000000)+tv.tv_usec;
}

void count_read(gzFile fq,uint32_t *i,const char *infile,FILE *out,double *mean_length){
	char *buf=(char *)calloc(1024,sizeof(char));
	uint32_t minLen=0,maxLen=0;
	uint64_t sum=0,sumQ1=0,sumQ2=0,sumFreq=0,*SeqLen=NULL,**Quality=NULL;
	Alloc2DArray(Quality,128,512);
	SeqLen=(uint64_t *)calloc(512,sizeof(uint64_t));
	while (gzgets(fq,buf,1024*sizeof(char))!=NULL){
		gzgets(fq,buf,1024*sizeof(char));
		uint16_t seqLen=strlen(buf)-1;
		SeqLen[seqLen]++;
		gzgets(fq,buf,1024*sizeof(char));
		gzgets(fq,buf,1024*sizeof(char));
		AssignQuality(Quality,buf,seqLen);
	}
	free(buf);
	statSeqLen(SeqLen,512,minLen,maxLen,*mean_length,sumFreq);
	*i=sumFreq;
	double Mean_Len=*mean_length/sumFreq;
	statQ(Quality,128,512,sum,53,sumQ1,63,sumQ2);
	Free2DArray(Quality,128);
	fprintf(out,"%s\t%lu\t%.0f\t%.0f\t%u\t%u\t%.3f\t%.3f\n",infile,sumFreq,*mean_length,Mean_Len,minLen,maxLen,1.0*sumQ1/sum*100,1.0*sumQ2/sum*100);
	if (globalArgs.LengthDetail){
		printSeqLen(SeqLen,512,out,minLen,maxLen,*mean_length);
	}
	free(SeqLen);
}

void display_usage(char * argv[]){
	char *buffer=(char* )malloc(10240*sizeof(char));
	const char* usage=
"\nCopyright (c) Berrygenomics 2014\n" \
"Contact: XiongXu <xuxiong@berrygenomics.com> <xiongxu@me.com> \n" \
"Discription:\n  This program is used for counting the read count with or without gzip compressed fastq files .\n" \
"Usage: %s file1.fq file2.fq ... [-o outfile] [-t thread] [-H] [-L] [-h] \n" \
"Example1:\n  %s /share/seq_dir1/Item/prenatal/140129_SN1116_0776_AC2WTGACXX/13C37198_L7_I012.R1.clean.fastq.gz \n" \
"\n" \
"   [-o OUTPUT] = OUTPUT file.default is stdout.                             [option]\n" \
"   [-H ]       = ouput the Header information.                              [option]\n" \
"   [-L ]       = ouput the read length in detail to the OUTPUT              [option]\n" \
"   [-t ]       = thread .default is the Infiles' number,if Infiles' number          \n" \
"                 > cpu cores,then the default is the cpu core number        [option]\n" \
"   [-h ] This helpful help screen.                                          [option]\n" \
"   Infiles     = Infile1,Infile2...                                         [required]\n" \
"\n";
	sprintf(buffer,usage,argv[0],argv[0]);
	fprintf(stderr,"%s",buffer);
	free(buffer);
	exit(1);
}

void load_fq(void *data,long i,int tid){
	fprintf(stderr,"tid: %d\n",tid);
	thread_arg_t *thread_arg=(thread_arg_t *)data;
	gzFile fq=open_input_stream(thread_arg->infiles[i]);
	count_read(fq,&thread_arg->read_count[i],thread_arg->infiles[i],thread_arg->out,&thread_arg->mean_len[i]);
	gzclose(fq);
}

int main(int argc, char *argv[]) {
	int opt = 0;
	globalArgs.infiles=NULL;
	globalArgs.numInfiles=0;
	globalArgs.header=0;
	globalArgs.LengthDetail=0;
	globalArgs.outfile="-";
	globalArgs.thread=CPU_NUM;
	const char *optString = "o:t:HLh?";
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
		switch( opt ) {
			case 'o':
				globalArgs.outfile = optarg;
				break;
			case 't':
				globalArgs.thread = atoi(optarg);
				break;
			case 'H':
				globalArgs.header++;
				break;
			case 'L':
				globalArgs.LengthDetail++;
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
	globalArgs.infiles = argv + optind;
	globalArgs.numInfiles = argc - optind;
	if (globalArgs.numInfiles<globalArgs.thread ) globalArgs.thread=globalArgs.numInfiles;

	long long begin=usec();
	
	FILE *out=fopen_output_stream(globalArgs.outfile);
	if (globalArgs.header) fprintf(out,"#Filename\tReadCount\tBaseCount\tMeanLen\tMinLen\tMaxLen\tQ20(%%)\tQ30(%%)\n");

	thread_arg_t *args=(thread_arg_t *)calloc(globalArgs.numInfiles,sizeof(thread_arg_t));
	args->infiles=globalArgs.infiles;
	args->mean_len=(double *)calloc(globalArgs.numInfiles,sizeof(double *));
	args->read_count=(uint32_t *)calloc(globalArgs.numInfiles,sizeof(uint32_t));
	args->out=out;

	kt_for(globalArgs.thread, load_fq, args, globalArgs.numInfiles);
	free(args->mean_len);
	free(args->read_count);
	free(args);

	fprintf(stderr,"Finished at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	fclose(out);
	return 0;
}
