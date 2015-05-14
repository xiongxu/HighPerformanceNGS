// gcc -g -O3 -Wall bam_count_GC.c  -o bam_count_GC -I./samtools-0.1.19 -L./samtools-0.1.19 -I/share/software/software/libgd_install/include -L/share/software/software/libgd_install/lib -lpthread -lgd -lpng -lz  -lfreetype -lm -lbam
#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <err.h>
#include <time.h>
#include <sys/time.h>
#include "sam.h"
#include "gd.h"
//#include "gdfontt.h"  //gdFontGetTiny()
//#include "gdfonts.h"//gdFontGetSmall()
//#include "gdfontmb.h" //gdFontGetMediumBold()
//#include "gdfontl.h"  //gdFontGetLarge()
#include "gdfontg.h"  //gdFontGetGiant()

typedef struct _data_buf_ {
	samfile_t *fp;
	FILE *out;
	uint32_t *windows;
	unsigned int **bins;
	unsigned int **len;
	float **GC;
	unsigned int total_mapped_reads;
	unsigned int *sum_count;
	unsigned long *sum_base;
	float *sum_GC;
	unsigned long n_count;
//	float **log_2_ratio;
}DataBuf;

struct globalArgs_t {
	char **infiles;
	unsigned short numInfiles;
	const char *outfile;
	int window;
	const char *region;
	int strand;
} globalArgs;

static inline float recursion(float bin[],int files_num);
void display_usage(char * argv[]);
static inline void cal_GC(const uint8_t *seq,int32_t l_qseq,unsigned short *n_GC);
static inline int fetch_func(const bam1_t *b, void *data);
void output_count_GC(DataBuf *databuf,const char *count_GC_prefix);
static inline uint32_t max(uint32_t *target_len,int32_t indices) ;
static inline float average(float data[], int n);
static inline float standard_deviation(float data[], int n);
static inline float recursion(float bin[],int files_num);
static inline void draw_hits(DataBuf *databuf,char *png_prefix);
static inline unsigned int get_quantile(unsigned int *array,unsigned int ubound,float quantile);
static inline int cmp(const void *a, const void *b);
static inline void calc_winGC(DataBuf *databuf);

static inline long long usec(void);


long long usec(void) {
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return (((long long)tv.tv_sec)*1000000)+tv.tv_usec;
}

void display_usage(char * argv[]){
	char *buffer=(char* )malloc(8092*sizeof(char));
	const char* usage=
"\nCopyright (c) 2012-2013\n" \
"Contact: XiongXu <xuxiong19880610@163.com> <xiongxu@me.com> \n" \
"Usage: %s [-o OUTFILE] [-w WINDOW_SIZE] [-r chr1:1-2000000] [-s 0] [-h] bamFile1 bamFile2 ..\n" \
"Discription:\n  This program is used for drawing the hits distribution in whole genome and calculating the reads count and GC content in each bin(specified by -w) with samtools indexed bam files.\n" \
"Example:\n  %s -o out -w 20000 /ifs5/ST_ANNO/USER/xuxiong/data/13X18312_L7_I021.R1.clean.fastq.sort.bam\n" \
"\n" \
"   [-o OUTPUT_FILE]  = OUTPUT file.                                   [required]\n" \
"   [-w WINDOW_SIZE]  = window size. default is 20000.                 [option]\n" \
"   [-r]              = region, default is whole genome.(chr1:1-20000) [option]\n" \
"   [-s]              = bool variant,strand or not, default is 0.      [option]\n" \
"   [-h]              = This helpful help screen.                      [option]\n" \
"   Infiles           = bam format input file(s),at least 1 bam file.  [required]\n" \
"\n";
	sprintf(buffer,usage,argv[0],argv[0]);
	fprintf(stderr,"%s",buffer);
	exit(1);
}

void cal_GC(const uint8_t *seq,int32_t l_qseq,unsigned short *n_GC) {
	int32_t i=0;
	int temp=0;
	for (i=0;i<l_qseq ;i++ ) {
		temp=bam1_seqi(seq, i);
		if (temp == 2 || temp == 4) (*n_GC)++;
	}
}

int fetch_func(const bam1_t *b, void *data) {
	DataBuf *databuf=(DataBuf *)data;
	const bam1_core_t *c = &b->core;
	if (b->core.tid < 0) return 0;
	if (c->flag&BAM_FUNMAP) return 0;
	if (*(databuf->windows+c->tid)==0) {
		*(databuf->windows+c->tid)=(unsigned int)databuf->fp->header->target_len[c->tid]/globalArgs.window+1;
		*(databuf->bins+c->tid)=(unsigned int *)calloc(*(databuf->windows+c->tid),sizeof(unsigned int));
		*(databuf->GC+c->tid)=(float *)calloc(*(databuf->windows+c->tid),sizeof(float));
		*(databuf->len+c->tid)=(unsigned int *)calloc(*(databuf->windows+c->tid),sizeof(unsigned int));
	}
	databuf->n_count++;
	uint32_t *cigar = bam1_cigar(b);
	int i, l;
	for (i = l = 0; i < c->n_cigar; ++i) {
		int op  = bam_cigar_op(cigar[i]);
		int len = bam_cigar_oplen(cigar[i]);
		if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
			l += len;
	}
//	fprintf(stdout,"%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%c\t%s\t%d\n",
//		databuf->fp->header->target_name[c->tid],c->pos, bam_calend(c,cigar) ,c->pos + l,l,bam_cigar2qlen(c,cigar), 
//		bam1_qname(b), c->qual, c->bin,c->mtid,c->mpos,c->isize,(c->flag&BAM_FREVERSE)? '-' : '+',bam_aux2Z(bam_aux_get(b,"MD")),bam_aux2i(bam_aux_get(b,"NH")));
	(databuf->total_mapped_reads)++;
	unsigned short current_window=(unsigned short)(c->pos/(globalArgs.window));
	unsigned short current_GC =0;
	cal_GC(bam1_seq(b),c->l_qseq,&current_GC);
	(*(*(databuf->bins+c->tid)+current_window))++;
	(*(*(databuf->GC+c->tid)+current_window))+=current_GC;
	(*(*(databuf->len+c->tid)+current_window))+=c->l_qseq;
	return 0;
}

void calc_winGC(DataBuf *databuf){
	int j;
	uint32_t k;
	for (j=0;j<databuf->fp->header->n_targets ;j++ ){
		for (k=0;k<*(databuf->windows+j) ;k++ ) {
			(*(databuf->sum_count+j))+=*(*(databuf->bins+j)+k);
			(*(databuf->sum_GC+j))+=*(*(databuf->GC+j)+k);
			(*(databuf->sum_base+j))+=*(*(databuf->len+j)+k);
			*(*(databuf->GC+j)+k)= *(*(databuf->GC+j)+k)!=0 ? *(*(databuf->GC+j)+k)/(*(*(databuf->len+j)+k))*100 : *(*(databuf->GC+j)+k);
		}
		*(databuf->sum_GC+j)=*(databuf->sum_GC+j)/(*(databuf->sum_base+j))*100;
	}
}

void output_count_GC(DataBuf *databuf,const char *count_GC_prefix){
	char *output=(char *)malloc(1024*sizeof(char));
	memset(output,0,1024*sizeof(char));
	sprintf(output,"%s.txt",count_GC_prefix);
	databuf->out=fopen(output,"wb");
	if(databuf->out == NULL) err(1, "Failed to open file (%s)", output);
	int j;
	uint32_t k;
	fprintf(databuf->out,"#chr\tchr_len\tchr_sum_read_count\tchr_sum_base\tchr_mean_cov\tchr_mean_GC%%");
	uint32_t max_target_window= max(databuf->fp->header->target_len,databuf->fp->header->n_targets)/(globalArgs.window)+1;
	for (k=0;k<max_target_window ;k++ ) {
		fprintf(databuf->out,"\t%u\tcount\tGC%%",k+1);
	}
	fprintf(databuf->out,"\n");
	for (j=0;j<databuf->fp->header->n_targets ;j++ ) {
		if (*(databuf->windows+j)==0) continue;
		fprintf(databuf->out,"%s\t%d\t%u\t%lu\t%f\t%f",databuf->fp->header->target_name[j],databuf->fp->header->target_len[j],*(databuf->sum_count+j),*(databuf->sum_base+j),(double)*(databuf->sum_base+j)/databuf->fp->header->target_len[j],*(databuf->sum_GC+j));
		for (k=0;k<*(databuf->windows+j) ;k++ ) {
			fprintf(databuf->out,"\t%d\t%u\t%f",k+1,*(*(databuf->bins+j)+k),*(*(databuf->GC+j)+k));
		}
		fprintf(databuf->out,"\n");
	}
	fclose(databuf->out);
	free(output);
}

/*
void output_count_GC(DataBuf *databuf,char *count_GC_prefix){
	int i,j,k;
	char *output=(char *)malloc(1024*sizeof(char));
	for (j=0;j<databuf->fp->header->n_targets ;j++ ) {
		if (*(databuf->windows+j)==0) continue;
		memset(output,0,1024*sizeof(char));
		sprintf(output,"%s_%s_RC_GC.txt",count_GC_prefix,databuf->fp->header->target_name[j]);
		FILE *out=fopen(output,"wb");
		if(out == NULL) err(1, "Failed to open file (%s)", output);
		fprintf(out,"%s\twindow",databuf->fp->header->target_name[j]);
		for (i=0;i<globalArgs.numInfiles ;i++ ) {
			fprintf(out,"\t%s_RC%%\t%s_GC%%",*(globalArgs.infiles+i),*(globalArgs.infiles+i));
		}
		fprintf(out,"\n#%s\t-",databuf->fp->header->target_name[j]);
		for (i=0;i<globalArgs.numInfiles ;i++ ) {
			fprintf(out,"\tsumRC(%u)\tmeanGC(%f)",*((databuf+i)->sum_count+j),*((databuf+i)->sum_GC+j));
		}
		fprintf(out,"\n");
		for (k=0;k<*(databuf->windows+j) ;k++ ) {
			fprintf(out,"%s\t%d",databuf->fp->header->target_name[j],k+1);
			for (i=0;i<globalArgs.numInfiles ;i++ ) {
				fprintf(out,"\t%u\t%f",*(*(databuf->bins+j)+k),*(*(databuf->GC+j)+k));
			}
			fprintf(out,"\n");
		}
		fclose(out);
	}
	free(output);
}
*/

uint32_t max(uint32_t *target_len,int32_t indices) {
	uint32_t max_len=0;
	unsigned short i=0;
	for (i=0;i<indices ;i++ ){
		max_len = *(target_len+i)>max_len ? *(target_len+i) : max_len;
	}
	return max_len;
}

float average(float data[], int n) {
	float mean=0.0;
	int i;
	for(i=0; i<n;++i){
		mean+=data[i];
	}
	return mean/n;
}

float standard_deviation(float data[], int n) {
	float sum_deviation=0.0;
	float mean=average(data,n);
	int i;
	for(i=0; i<n;++i){
		sum_deviation+=(data[i]-mean)*(data[i]-mean);
	}
	return sqrt(sum_deviation/(n-1));
}

float recursion(float bin[],int files_num) {
//	qsort((void *)bin,files_num,sizeof(int),comp);
	float sd=standard_deviation(bin,files_num);
	float mean= average(bin,files_num);
	float upper=mean+3*sd;
	float down=mean-3*sd;
	down=down>0?down:0;
	
	int flag=0;
	int i=0,j=0;
	for (i=0;i<files_num ;i++ ){
		if (bin[i]<down||bin[i] > upper){
			flag=1;
//			fprintf(stdout,"%d %f\n",i,bin[i]);
			for (j=i;j<files_num-1 ;j++){
				bin[j]=bin[j+1];
			}
			files_num--;
		}
	}
	if (flag) {
//		fprintf(stdout,"SD: %f\nMean: %f\nMean-theta*SD: %f\nMean+theta*SD: %f\n\n",sd,mean,down,upper);
//		print_array(bin,files_num,stdout);
		recursion(bin,files_num);
	}
	else{
//		print_array(bin,files_num,stdout);
//		fprintf(stdout,"previous mean: %f\n",mean);
		mean= average(bin,files_num);
//		fprintf(stdout,"SD: %f\nLast mean: %f\nMean-theta*SD: %f\nMean+theta*SD: %f\n",sd,mean,down,upper);
		return mean;
	}
}

int cmp(const void *a, const void *b){
	return *(int *)a-*(int *)b;
}

unsigned int get_quantile(unsigned int *array,unsigned int ubound,float quantile){
	unsigned int *temp_array=calloc(ubound,sizeof(unsigned int));
	memcpy(temp_array,array,ubound*sizeof(unsigned int));
	qsort(temp_array,ubound,sizeof(unsigned int),cmp);
	unsigned int quantile_value=*(temp_array+(int)(quantile*ubound));
//	fprintf(stderr,"%d\t%d\t%d\n",ubound,quantile_value,*(temp_array+ubound-1));
	free(temp_array);
	return quantile_value;
}

void draw_hits(DataBuf *databuf,char *png_prefix){
	int i=0,j=0;
	uint32_t k=0;
	char *output=(char *)malloc(1024*sizeof(char));
	memset(output,0,1024*sizeof(char));
	sprintf(output,"%s_hits.png",png_prefix);
	databuf->out=fopen(output,"wb");
	if(databuf->out == NULL) err(1, "Failed to open file (%s)", output);

	uint32_t max_target_window= max(databuf->fp->header->target_len,databuf->fp->header->n_targets)/(globalArgs.window)+1;
	unsigned int left=50;
	unsigned int top=50;
	float each_window_size = 0.2;
	unsigned int each_chromosome_height = 40;
	unsigned int Height = databuf->fp->header->n_targets * (each_chromosome_height+10) + 2*top;
	unsigned int Width = max_target_window*each_window_size + 2*left;
	
	gdImagePtr im = gdImageCreate(Width,Height+top);
	int color[12];
	int rgb[12][3]={
			{255,255,255},		/*******set the background_color*******/
			{70, 130, 180},{255, 140, 0},{160, 82, 45},{135, 206, 235},{107, 142, 35},{106, 90, 205},
			{119, 136, 153},{218, 165, 32},{178, 34, 34},{255, 0, 255},{0, 255, 255}
	};
	for (i=0;i<12 ;i++ ) {
		color[i] = gdImageColorExact(im,rgb[i][0], rgb[i][1],rgb[i][2]);
		if (color[i] != (-1)) gdImageColorDeallocate(im, color[i]);
		color[i] = gdImageColorAllocate(im,rgb[i][0], rgb[i][1],rgb[i][2]);
	}
	int im_black = gdImageColorAllocate(im,0,0,0);
//	gdImageRectangle(im, left,top,Width-left,Height-top, im_black);

	float x=left;
	unsigned int y=top;
	for (j=0;j<databuf->fp->header->n_targets ;j++ ) {
//		*(databuf->log_2_ratio+j)=(float *)calloc(*(databuf->windows+j),sizeof(float));
		x=left;
		y+=each_chromosome_height+10;
		gdImageString(im,gdFontGetGiant(),0,y-20,(unsigned char *)databuf->fp->header->target_name[j],im_black);
		gdImageRectangle(im, (unsigned int)x,y-each_chromosome_height,left+*(databuf->windows+j)*each_window_size,y, im_black);
		if (*(databuf->windows+j)==0) continue;
		unsigned int percent95_value= get_quantile(*(databuf->bins+j),*(databuf->windows+j),0.95);
		for (k=0;k<*(databuf->windows+j) ;k++ ) {
			x+=each_window_size;
			gdImageSetPixel(im,(unsigned int) x, y-(percent95_value==0?0:*(*(databuf->bins+j)+k)*each_chromosome_height/percent95_value), color[1]);
		}
	}
	x=left;
	gdImageFilledRectangle(im, (unsigned int)x,y,(unsigned int)x+each_chromosome_height,y+each_chromosome_height, color[1]);
	gdImageString(im,gdFontGetGiant(),(unsigned int)x+left,y+20,(unsigned char *)png_prefix,color[1]);
	y+=each_chromosome_height+10;
	gdImagePng(im,databuf->out);
	gdImageDestroy(im);
	fclose(databuf->out);
	free(output);
}

int main(int argc, char *argv[])
{
	int opt = 0;
	globalArgs.infiles=NULL;
	globalArgs.numInfiles=0;
	globalArgs.outfile="out";
	globalArgs.window=20000;
	globalArgs.region="-";
	globalArgs.strand=0;
	const char *optString = "o:w:r:s:h?";
	if (argc<2) display_usage(argv);
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
		switch( opt ) {
			case 'o':
				globalArgs.outfile = optarg;
				break;
			case 'w':
				globalArgs.window = atoi(optarg);
				break;
			case 'r':
				globalArgs.region = optarg;
			case 's':
				globalArgs.strand = atoi(optarg);
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

	DataBuf *databuf=(DataBuf *)calloc(globalArgs.numInfiles,sizeof(DataBuf));

	uint32_t i=0;
	long long begin;
	begin=usec();
	for (i=0;i<globalArgs.numInfiles ;i++ ) {
		if (((databuf+i)->fp = samopen(*(globalArgs.infiles+i), "rb", 0)) == 0) {
			err(1, "bam2bed: Fail to open BAM file %s\n", *(globalArgs.infiles+i));
		}
		(databuf+i)->total_mapped_reads=0;
		(databuf+i)->bins = (unsigned int **)calloc((databuf+i)->fp->header->n_targets,sizeof(unsigned int *));
		(databuf+i)->GC = (float **)calloc((databuf+i)->fp->header->n_targets,sizeof(float *));
		(databuf+i)->len = (unsigned int **)calloc((databuf+i)->fp->header->n_targets,sizeof(unsigned int *));
//		(databuf+i)->log_2_ratio = (float **)calloc((databuf+i)->fp->header->n_targets,sizeof(float *));
		(databuf+i)->windows = (uint32_t *)calloc((databuf+i)->fp->header->n_targets,sizeof(uint32_t));
		(databuf+i)->sum_count = (unsigned int *)calloc((databuf+i)->fp->header->n_targets,sizeof(unsigned int));
		(databuf+i)->sum_GC = (float *)calloc((databuf+i)->fp->header->n_targets,sizeof(float));
		(databuf+i)->sum_base = (unsigned long *)calloc((databuf+i)->fp->header->n_targets,sizeof(unsigned long));
//		uint32_t max_target_window= max((databuf+i)->fp->header->target_len,(databuf+i)->fp->header->n_targets)/globalArgs.window+1;
//		fprintf(stderr,"%u|%d\n",max_target_window,(databuf+i)->fp->header->n_targets);
		if (strncmp(globalArgs.region,"-",6)==0) { /* if a region is not specified */
			bam1_t *b = bam_init1();
			while (samread((databuf+i)->fp, b) >= 0) fetch_func(b, databuf+i);
			bam_destroy1(b);
		}
		else {
			int ref, beg, end;
			bam_index_t *idx;
			if ((idx = bam_index_load(*(globalArgs.infiles+i))) == 0) {
				fprintf(stderr, "bam2bed: BAM indexing file is not available.\n");
				return 1;
			}
			bam_parse_region((databuf+i)->fp->header, globalArgs.region, &ref, &beg, &end);
			if (ref < 0) {
				fprintf(stderr, "bam2bed: Invalid region %s\n", globalArgs.region);
				return 1;
			}
			fprintf(stdout,"%s\t%d\t%d\n",(databuf+i)->fp->header->target_name[ref],beg,end);
			bam_fetch((databuf+i)->fp->x.bam, idx, ref, beg, end, (databuf+i), fetch_func);
			bam_index_destroy(idx);
		}
		fprintf(stderr,"Done load bam file %s at %.3f s\n",*(globalArgs.infiles+i),(double)(usec()-begin)/CLOCKS_PER_SEC);
		calc_winGC(databuf+i);
		draw_hits(databuf+i,*(globalArgs.infiles+i));
		fprintf(stderr,"Done draw hit %s_hits.png at %.3f s\n",*(globalArgs.infiles+i),(double)(usec()-begin)/CLOCKS_PER_SEC);
		fprintf(stderr,"%lu\n",databuf->n_count);
	}
	output_count_GC(databuf,globalArgs.outfile);
	fprintf(stderr,"Done output %s.txt at %.3f s\n",globalArgs.outfile,(double)(usec()-begin)/CLOCKS_PER_SEC);
	
	int j=0;
	for (i=0;i<globalArgs.numInfiles ;i++ ) {
		for (j=0;j<(databuf+i)->fp->header->n_targets ;j++ ) {
			free(*((databuf+i)->bins+j));
			free(*((databuf+i)->GC+j));
		}
		free((databuf+i)->windows);
		free((databuf+i)->sum_count);
		free((databuf+i)->sum_GC);
		free((databuf+i)->sum_base);
		free((databuf+i)->bins);
		free((databuf+i)->GC);
		samclose((databuf+i)->fp);
	}
	free(databuf);
	return 0;
}
