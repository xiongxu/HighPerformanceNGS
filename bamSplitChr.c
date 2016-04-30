// gcc -g -O3 -Wall  bamSplit.c -o bamSplit -I./samtools-0.1.19 -I./zlib-1.2.8/include -L./samtools-0.1.19 -L./zlib-1.2.8/lib -lpthread -lm -lz -lbam
#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <err.h>
#include <time.h>
#include <sys/time.h>
#include <libgen.h>
#include "sam.h"

typedef struct _data_buf_ {
	samfile_t *fp;
	samfile_t *out;
	uint32_t count;
}DataBuf;

struct globalArgs_t {
	char **infiles;
	unsigned short numInfiles;
	const char *outfile;	
} globalArgs;

static inline long long usec(void);
static inline void bam_fetch_chr(bam_index_t *idx,const char *region,DataBuf *databuf);
static inline bam_index_t *load_bam_index(char *filename);
static inline int view_func(const bam1_t *b, void *data);
void display_usage(char * argv[]);

long long usec(void) {
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return (((long long)tv.tv_sec)*1000000)+tv.tv_usec;
}

void display_usage(char * argv[]){
	char *buffer=(char* )malloc(8092*sizeof(char));
	const char* usage=
"\nCopyright (c) 2016\n" \
"Contact: XiongXu <xuxiong19880610@163.com> <xiongxu@me.com> \n" \
"Usage: %s [-o OUTFILE_PREFIX] [-h] bamFile1 bamFile2 ..\n" \
"Discription:\n  This program is used for spliting bam files into some chromosomes'bam.\n" \
"Example1:\n  %s /share/work1/staff/xuxiong/lvliyaELL_L1_I008_tophat_mis2/accepted_hits.bam\n" \
"\n" \
"   [-o OUTPUT_FILE]  = OUTPUT file.                                   [required]\n" \
"   -1                   use fast BAM compression                      [option]\n" \
"   -u                   uncompressed BAM output                       [option]\n" \
"   [-h]              = This helpful help screen.                      [option]\n" \
"   Infiles           = bam format input file(s),at least 1 bam file.  [required]\n" \
"\n";
	sprintf(buffer,usage,argv[0],argv[0]);
	fprintf(stderr,"%s",buffer);
	exit(1);
}

int view_func(const bam1_t *b, void *data) {
	DataBuf *databuf = (DataBuf *)data;
	++databuf->count;
	// if (!process_aln(((DataBuf *)data)->out->header, (bam1_t*)b))
	samwrite(databuf->out, b);
	return 0;
}

bam_index_t *load_bam_index(char *filename) {
	bam_index_t *idx;
	if ((idx = bam_index_load(filename)) == 0) {
		fprintf(stderr, "bam2bed: BAM indexing file is not available.\n");
		exit(1);
	}
	return idx;
}

void bam_fetch_chr(bam_index_t *idx,const char *region,DataBuf *databuf) {
	int ref, beg, end;
	bam_parse_region(databuf->fp->header,region, &ref, &beg, &end);
	if (ref < 0) {
		fprintf(stderr, "Invalid region %s\n",region);
		exit(1);
	}
	bam_fetch(databuf->fp->x.bam, idx, ref, beg, end, databuf, view_func);
}

int main(int argc, char *argv[]) {
	int opt = 0;
	globalArgs.infiles=NULL;
	globalArgs.numInfiles=0;
	globalArgs.outfile=NULL;
	int compress_level=-1;
	const char *optString = "o:w:r:s:u:1:h?";
	if (argc<2) display_usage(argv);
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
		switch( opt ) {
			case 'o':
				globalArgs.outfile = optarg;
				break;
			case 'u': compress_level = 0; break;
			case '1': compress_level = 1; break;
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
	
	long long begin=usec();

	DataBuf *databuf=(DataBuf *)calloc(globalArgs.numInfiles,sizeof(DataBuf));
	uint32_t i=0;
	int j=0;
	char out_mode[5], fn_out[256];
	memcpy(out_mode,"wbh",3);
	if (compress_level >= 0) {
		char tmp[2];
		tmp[0] = compress_level + '0'; tmp[1] = '\0';
		strcat(out_mode, tmp);
	}
	for (i=0;i<globalArgs.numInfiles ;i++ ) {
		if (!globalArgs.outfile) globalArgs.outfile = globalArgs.infiles[i];
		if (((databuf+i)->fp = samopen(*(globalArgs.infiles+i), "rb", 0)) == 0) {
			err(1, "bam2bed: Fail to open BAM file %s\n", *(globalArgs.infiles+i));
		}
		bam_index_t *idx=load_bam_index(*(globalArgs.infiles+i));
		for (j=0;j<(databuf+i)->fp->header->n_targets ;j++ ) {
			sprintf(fn_out,"%s_%s.bam",globalArgs.outfile,(databuf+i)->fp->header->target_name[j]);
			if ( ((databuf+i)->out = samopen(fn_out, out_mode, (databuf+i)->fp->header)) == 0 ){
				err(1, "[main_samview] fail to open \"%s\" for writing.\n", fn_out);
			}
			(databuf+i)->count=0;
			bam_fetch_chr(idx,(databuf+i)->fp->header->target_name[j],databuf+i);
			fprintf(stderr,"chr: %s\tchr_len: %d\treads_count: %u at %.3f s\n",(databuf+i)->fp->header->target_name[j],(databuf+i)->fp->header->target_len[j],(databuf+i)->count,(double)(usec()-begin)/CLOCKS_PER_SEC);
			samclose((databuf+i)->out);
		}
		bam_index_destroy(idx);
		samclose((databuf+i)->fp);
		globalArgs.outfile =NULL;
		fprintf(stderr,"splited %s into each chromosome at %.3f s\n",*(globalArgs.infiles+i),(double)(usec()-begin)/CLOCKS_PER_SEC);
	}
	free(databuf);
	return 0;
}
