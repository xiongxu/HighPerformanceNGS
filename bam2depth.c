// gcc -g -O3 -Wall bam2depth.c -o bam2depth -I./samtools-0.1.19 -I./zlib-1.2.8/include -L./samtools-0.1.19 -L./zlib-1.2.8/lib -I. -L. -lhashtbl -lpthread -lm -lz -lbam
#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <err.h>
#include <time.h>
#include <sys/time.h>
#include <libgen.h>
#include "sam.h"
#include "hashtbl.h"
#include "IO_stream.h"

#define HASH_SIZE 5000000

typedef struct _data_buf_ {
	samfile_t *fp;
	HASHTBL *Start;
	HASHTBL *End;
}DataBuf;

struct globalArgs_t {
	char **infiles;
	unsigned short numInfiles;
	const char *outfile;
	uint32_t window;
	const char *region;
	int strand;
	int wig;
} globalArgs;

static inline long long usec(void);
void output_bins(DataBuf *databuf,int ref,double *bins,int windows,FILE *stream);
void hash2BedGraph(DataBuf *databuf,int ref,double *bins,int windows,FILE *bedGraph);
static inline void hash_data_delete(void *DATA);
static inline int chr2int(const char *str);
static inline void overlap(int *j,int ref,int *subject_count,DataBuf *databuf,int windows,double *bins,uint32_t last_start,uint32_t last_end,double last_depth);
void bam_fetch_chr(bam_index_t *idx,const char *region,DataBuf *databuf);
bam_index_t *load_bam_index(char *filename);
static inline int fetch_func(const bam1_t *b, void *data);
static inline int insertInt2Hash(HASHTBL *Start,int pos,char *value);
void display_usage(char * argv[]);

long long usec(void) {
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return (((long long)tv.tv_sec)*1000000)+tv.tv_usec;
}

void display_usage(char * argv[]){
	char *buffer=(char* )malloc(8092*sizeof(char));
	const char* usage=
"\nCopyright (c) 2015\n" \
"Contact: XiongXu <xuxiong19880610@163.com> <xiongxu@me.com> \n" \
"Usage: %s [-o OUTFILE] [-w WINDOW_SIZE] [-g gff_file] [-r chr1:1-2000000] [-W is wig] [-s 0] [-h] bamFile1 bamFile2 ..\n" \
"Discription:\n  This program is used for converting bam files to bedgraph format files,meanwhile calculate each window's base count and mean depth.\n" \
"Example1:\n  %s -o out -w 20000 /share/work1/staff/xuxiong/lvliyaELL_L1_I008_tophat_mis2/accepted_hits.bam > depth \n" \
"Example2:\n  %s -o out -w 20000 /share/work1/staff/xuxiong/twins/batch2/test_20131224/RT144FD_L1_I001.R1.clean.fastq.gz_1224/accepted_hits.bam > depth\n"\
"\n" \
"   [-o OUTPUT_FILE]  = OUTPUT file.                                   [required]\n" \
"   [-w WINDOW_SIZE]  = window size. default is 20000.                 [option]\n" \
"   [-W is wig]  = whether to generate wig format file. default not.   [option]\n" \
"   [-r]              = region, default is whole genome.(chr1:1-20000) [option]\n" \
"   [-s]              = bool variant,strand or not, default is 0.      [option]\n" \
"   [-h]              = This helpful help screen.                      [option]\n" \
"   Infiles           = bam format input file(s),at least 1 bam file.  [required]\n" \
"\n";
	sprintf(buffer,usage,argv[0],argv[0],argv[0]);
	fprintf(stderr,"%s",buffer);
	exit(1);
}

int insertInt2Hash(HASHTBL *Start,int pos,char *value) {
	int2char(pos,(UBYTE *)value);
	void *val_data=hashtbl_get(Start,value);
	if (val_data==NULL) {
		int *init=(int *)malloc(sizeof(int));
		*init=1;
		hashtbl_insert(Start, value, (void *)init);
	}
	else{
		(*(int *)val_data)++;
	}
	return 0;
}

int fetch_func(const bam1_t *b, void *data) {
	DataBuf *databuf=(DataBuf *)data;
	uint32_t *cigar = bam1_cigar(b);
	const bam1_core_t *c = &b->core;
	if (c->flag&BAM_DEF_MASK || b->core.tid < 0) return 0;
	int i, l;
	char *buf=(char *)malloc(32*sizeof(char));
	unsigned int temp_start= c->pos;
	for (i = l = 0; i < c->n_cigar; ++i) {
		int op = bam_cigar_op(cigar[i]);
		if (op == BAM_CINS) continue;
		//if ( op == BAM_CDEL || op == BAM_CREF_SKIP || op ==BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
		if ( op == BAM_CDEL || op == BAM_CREF_SKIP) {
			l = bam_cigar_oplen(cigar[i]);
			temp_start+=l;
		}else if (op == BAM_CMATCH ) {
			l = bam_cigar_oplen(cigar[i]);
			insertInt2Hash(databuf->Start,temp_start,buf);
			temp_start+=l;
			insertInt2Hash(databuf->End,temp_start,buf);
		}
	}
	free(buf);
	return 0;
}

bam_index_t *load_bam_index(char *filename){
	bam_index_t *idx;
	if ((idx = bam_index_load(filename)) == 0) {
		fprintf(stderr, "bam2bed: BAM indexing file is not available.\n");
		exit(1);
	}
	return idx;
}

void bam_fetch_chr(bam_index_t *idx,const char *region,DataBuf *databuf){
	int ref, beg, end;
	bam_parse_region(databuf->fp->header,region, &ref, &beg, &end);
	if (ref < 0) {
		fprintf(stderr, "bam2bed: Invalid region %s\n", globalArgs.region);
		exit(1);
	}
//	fprintf(stdout,"%s\t%d\t%d\n",databuf->fp->header->target_name[ref],beg,end);
	bam_fetch(databuf->fp->x.bam, idx, ref, beg, end, databuf, fetch_func);
}

void overlap(int *j,int ref,int *subject_count,DataBuf *databuf,int windows,double *bins,uint32_t last_start,uint32_t last_end,double last_depth){
	if (*subject_count>1) {//deal whith the overlapping query
		if((*j-*subject_count)>=0){
			(*j)-=*subject_count ;
		}else{
			*j=0;
		}
	}
	uint32_t window_start,window_end;
	*subject_count=0;
	while (*j<=windows){
		window_start=globalArgs.window*(*j);
		window_end= (*j+1)*globalArgs.window;
		if (window_end>databuf->fp->header->target_len[ref]){
			window_end=databuf->fp->header->target_len[ref];
		}
		if (last_end<window_start) {
			break;
		}else{
			if (last_start<=window_start) {
				if (last_end<=window_end) {
					bins[*j]+=(last_end-window_start)*last_depth;
					(*subject_count)++;
					break;
				}else{
					bins[(*j)++]+=(window_end-window_start)*last_depth;
					(*subject_count)++;
				}
			}else{
				if (last_start<window_end) {
					if (last_end<=window_end) {
						bins[*j]+=(last_end-last_start)*last_depth;
						(*subject_count)++;
						break;
					}else{
						bins[(*j)++]+=(window_end-last_start)*last_depth;
						(*subject_count)++;
					}
				}else{
					(*j)++;
				}
			}
		}
	}
}

int chr2int(const char *str){
	int temp = 0;
	const char *ptr = str;
	if (*str == '-' || *str == '+') str++;
	while(*str){
		if (*str == 'X') {
			return 23;
		}else if (*str == 'Y') {
			return 24;
		}else if (*str == 'M') {
			return 25;
		}
		if ((*str < '0') || (*str > '9')) {str++;continue;}
		temp = temp * 10 + (*str++ - '0');
	}
	if (*ptr == '-') temp = -temp;
	return temp;
}

void hash_data_delete(void *DATA) {
	ENTRY *node=(ENTRY *)DATA;
	free(node->data);
	free(node);
}

void hash2BedGraph(DataBuf *databuf,int ref,double *bins,int windows,FILE *bedGraph){
	int all_keys_count=databuf->Start->count+databuf->End->count;
	if (!all_keys_count) return;
	int i=0,j=0,prevkey=0,last_start=0,last_end=0,last_depth=0,Count=0,subject_count=0,pos=0;
	char **all_keys=union_hashed_keys(databuf->Start,databuf->End);
	for (i = 0; i < all_keys_count; ++i) {
		char2int((UBYTE *)all_keys[i],&pos);
		if (pos && pos==prevkey) continue;
		if (prevkey) {
			if (last_depth==Count) {
				prevkey=last_start;
			}
			else{
				if (last_depth){
					fprintf(bedGraph,"%s\t%d\t%d\t%d\n",databuf->fp->header->target_name[ref],last_start,last_end,last_depth);
					overlap(&j,ref,&subject_count,databuf,windows,bins,last_start,last_end,(double)last_depth);
				}
			}
		}
		last_start=prevkey;
		last_end=pos;
		last_depth=Count;
		void *data=hashtbl_get(databuf->Start,all_keys[i]);
		if (data!=NULL) Count+= *(int*)data;
		data=hashtbl_get(databuf->End,all_keys[i]);
		if (data!=NULL) Count-= *(int*)data;
		prevkey=pos;
	}
	if (last_depth){
		overlap(&j,ref,&subject_count,databuf,windows,bins,last_start,last_end,(double)last_depth);
		fprintf(bedGraph,"%s\t%d\t%d\t%d\n",databuf->fp->header->target_name[ref],last_start,last_end,last_depth);
	}
	free(all_keys);
}

void output_bins(DataBuf *databuf,int ref,double *bins,int windows,FILE *stream){
	int i;
	for (i=0;i<windows ;i++ ){
		int window_start = globalArgs.window*i;
		int window_end = globalArgs.window*(i+1) >databuf->fp->header->target_len[ref] ? databuf->fp->header->target_len[ref] : globalArgs.window*(i+1);
		fprintf(stream,"%s\t%d\t%d\t%.2f\n",
			databuf->fp->header->target_name[ref],window_start,window_end,bins[i]/globalArgs.window);
	}
}

void output_bins_wig(DataBuf *databuf,int ref,double *bins,int windows,FILE *stream){
	int i;
	fprintf(stream,"variableStep chrom=%s span=%d\n",databuf->fp->header->target_name[ref],globalArgs.window);
	for (i=0;i<windows ;i++ ){
		int window_start = globalArgs.window*i;
		if (bins[i]) fprintf(stream,"%d\t%.2f\n",window_start,bins[i]/globalArgs.window);
	}
}

int main(int argc, char *argv[])
{
	int opt = 0;
	globalArgs.infiles=NULL;
	globalArgs.numInfiles=0;
	globalArgs.outfile="-";
	globalArgs.window=20000;
	globalArgs.region="-";
	globalArgs.strand=0;
	globalArgs.wig=0;
	const char *optString = "o:w:r:s:Wh?";
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
			case 'W':
				globalArgs.wig++;
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
	
	long long begin;
	begin=usec();

	DataBuf *databuf=(DataBuf *)calloc(globalArgs.numInfiles,sizeof(DataBuf));
	uint32_t i=0;
	int j=0;
	char *suffix=(char *)calloc(64,sizeof(char));

	for (i=0;i<globalArgs.numInfiles ;i++ ) {
		if (((databuf+i)->fp = samopen(*(globalArgs.infiles+i), "rb", 0)) == 0) {
			err(1, "bam2bed: Fail to open BAM file %s\n", *(globalArgs.infiles+i));
		}
		FILE *WIG=NULL,*bedGraph=NULL,*depth=NULL,*chrSize=NULL;
		sprintf(suffix,".%u.bedGraph",i+1);
		bedGraph=fcreat_outfile(basename(*(globalArgs.infiles+i)),suffix);
		sprintf(suffix,".%u.depth",i+1);
		depth=fcreat_outfile(globalArgs.outfile,suffix);
		if (globalArgs.wig){
			sprintf(suffix,".%u.wig",i+1);
			WIG=fcreat_outfile(globalArgs.outfile,suffix);
			sprintf(suffix,".%u.chromSize.txt",i+1);
			chrSize=fcreat_outfile(globalArgs.outfile,suffix);
		}
		bam_index_t *idx=load_bam_index(*(globalArgs.infiles+i));
		double **Bins=(double **)calloc((databuf+i)->fp->header->n_targets,sizeof(double *));
		int *Windows=(int *)calloc((databuf+i)->fp->header->n_targets,sizeof(int));
		for (j=0;j<(databuf+i)->fp->header->n_targets ;j++ ) {
			Windows[j]=(databuf+i)->fp->header->target_len[j]/globalArgs.window+1;
			(databuf+i)->Start=hashtbl_create(HASH_SIZE, NULL);
			(databuf+i)->End=hashtbl_create(HASH_SIZE, NULL);
			bam_fetch_chr(idx,(databuf+i)->fp->header->target_name[j],databuf+i);
			Bins[j]=(double *)calloc(Windows[j],sizeof(double));
			hash2BedGraph(databuf+i,j,Bins[j],Windows[j],bedGraph);
			hashtbl_destroy((databuf+i)->Start,hash_data_delete);
			hashtbl_destroy((databuf+i)->End,hash_data_delete);
			output_bins(databuf+i,j,Bins[j],Windows[j],depth);
			if (globalArgs.wig) output_bins_wig(databuf+i,j,Bins[j],Windows[j],WIG);
			free(Bins[j]);
			if (globalArgs.wig) fprintf(chrSize,"%s\t%d\n",(databuf+i)->fp->header->target_name[j],(databuf+i)->fp->header->target_len[j]);
			fprintf(stderr,"%s at %.3f s\n",(databuf+i)->fp->header->target_name[j],(double)(usec()-begin)/CLOCKS_PER_SEC);
		}
		free(Windows);
		free(Bins);
		bam_index_destroy(idx);
		samclose((databuf+i)->fp);
		fclose(bedGraph);
		fclose(depth);
		if (globalArgs.wig){
			fclose(WIG);
			fclose(chrSize);
		}
		fprintf(stderr,"Converted %s to wig format at %.3f s\n",*(globalArgs.infiles+i),(double)(usec()-begin)/CLOCKS_PER_SEC);
	}
	free(suffix);
	free(databuf);
	return 0;
}
