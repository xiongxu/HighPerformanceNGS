//gcc -shared -std=gnu99 -fPIC -g -O2 Rgzfastq_uniq_3.c -o Rgzfastq_uniq_3.so -I/share/software/software/R-3.0_install/R-3.0.1/include -DNDEBUG  -I/usr/local/include -I/share/software/software/zlib_1.2.8_install/include -I/home/xuxiong/work/c -L/usr/local/lib64 -L/share/software/software/R-3.0_install/R-3.0.1/lib -L/home/xuxiong/work/c -L/share/software/software/zlib_1.2.8_install/lib -lR -lz -lhashtbl
/*
dyn.load("/home/xuxiong/work/R/Rgzfastq_uniq.so");
List<-.Call("qsort_hash_count","/share/work1/staff/xuxiong/test/13C37198_L7_I012.R1.clean.fastq.gz","")
M<-List[[3]]
M<-M[apply(M,1,function(X){!all(X==0)}),]
matrix.axes <- function(data) {
	x <- (1:dim(data)[1] - 1) / (dim(data)[1] - 1);
	axis(side=1, at=x, labels=1:dim(data)[1], las=2);
	x <- (1:dim(data)[2] - 1) / (dim(data)[2] - 1);
	axis(side=2, at=x, labels=1:dim(data)[2], las=2);
	grid(nx=(dim(data)[1]-1), ny=(dim(data)[2]-1), col="black", lty="solid");
}
filled.contour(M, plot.axes=matrix.axes(M), col=colorpanel(50, "white", "red"), nlevels=50, main="Protein-Protein Interaction Potential")
*/

#include <stdio.h>
#include <R.h>
#include <Rdefines.h>
#include <zlib.h>
#include <time.h>
#include <sys/time.h>
#include "hashtbl.h"

#define ELECNT 10000000
#define MaxLen 300
#define free_Fastq(line)	\
	do{							\
		free((line)->seq);		\
		free((line)->quality);	\
		free((line));			\
	}while (0)

#define Alloc2DArray(Array,nrow,ncol,Array1,nrow1,ncol1,Array2,Array2Count,Array3,Array3Count)	\
do{															\
	Array=(int *)calloc((nrow)*(ncol),sizeof(int));			\
	Array1=(int *)calloc((nrow1)*(ncol1),sizeof(int));		\
	Array2=(int *)calloc(Array2Count,sizeof(int));			\
	Array3=(double *)calloc(Array3Count,sizeof(double));	\
}while (0)

#define AssignQuality(Array,Q)	\
do{												\
	int index=0;								\
	for(;*(Q+index);index++){					\
		(*(Array+(int)*(Q+index)+128*index))++;	\
	}											\
}while (0)

#define STATSEQ(seq,GC,Nuc,L)	\
do{												\
	for(GC=0,L=0;*(seq+L);L++) {				\
		if (seq[L]=='G' || seq[L]=='C') (GC)++;	\
		(*(Nuc+5*L+ntVal[seq[L]]))++;			\
	}											\
	(GC)/=L;									\
}while (0)

#define LOAD_HASHTBL(fq1,fq2)	\
	initNtVal();									\
	HASHTBL *hashtbl=NULL;							\
	const char *fq_file1 = CHARACTER_VALUE(fq1);	\
	const char *fq_file2 = CHARACTER_VALUE(fq2);	\
	unsigned long total_reads_num=0;				\
	long long begin=usec();							\
	int **Quality=(int **)calloc((strcmp(fq_file2,"")?2:1),sizeof(int *));	\
	int **Length=(int **)calloc((strcmp(fq_file2,"")?2:1),sizeof(int *));	\
	int **Nucleotide=(int **)calloc((strcmp(fq_file2,"")?2:1),sizeof(int *));	\
	double **PE_GC = load_fastq_file(fq_file1,fq_file2,&total_reads_num,begin,&hashtbl,Quality,Length,Nucleotide)

typedef struct _fastq_
{
//	char *name;
	char *seq;
	char *quality;
}Fastq;

static inline Fastq *readNextNode(gzFile fq,char *buf);
gzFile creat_gzinfile(const char *infile);
double **load_fastq_file(const char *fq_file1,const char *fq_file2,unsigned long *total_reads_count,long long begin,HASHTBL **tbl,int **Quality,int **Length,int **Nucleotide);
static inline int compare_hashed_data_count(const void *a, const void *b);
SEXP qsort_hash_count(SEXP fq1,SEXP fq2);
static inline long long usec(void);
static inline void int2str(int c, int base, char *ret);
void initNtVal(void);

int ntVal[256];
//char valToNt[5];
//char revNt[256];
#define T_BASE_VAL 0
#define U_BASE_VAL 0
#define C_BASE_VAL 1
#define A_BASE_VAL 2
#define G_BASE_VAL 3
#define N_BASE_VAL 4   /* Used in 1/2 byte representation. */

void initNtVal(void){
	int i;
	for (i=0; i<256; i++){
		ntVal[i] = T_BASE_VAL;
//		revNt[i]=0;
	}
	ntVal['t'] = ntVal['T'] = T_BASE_VAL;
	ntVal['u'] = ntVal['U'] = U_BASE_VAL;
	ntVal['c'] = ntVal['C'] = C_BASE_VAL;
	ntVal['a'] = ntVal['A'] = A_BASE_VAL;
	ntVal['g'] = ntVal['G'] = G_BASE_VAL;
	ntVal['.'] = ntVal['N'] = N_BASE_VAL;
//	valToNt[T_BASE_VAL]= revNt['A']=revNt['a'] = 'T';
//	valToNt[C_BASE_VAL]= revNt['G']=revNt['g'] = 'C';
//	valToNt[A_BASE_VAL]= revNt['T']=revNt['t'] = 'A';
//	valToNt[G_BASE_VAL]= revNt['C']=revNt['c'] = 'G';
//	valToNt[N_BASE_VAL]= revNt['N']=revNt['.'] = 'N';
}

long long usec(void){
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return (((long long)tv.tv_sec)*1000000)+tv.tv_usec;
}

Fastq *readNextNode(gzFile fq,char *buf) {
	Fastq *line = (Fastq *)malloc(sizeof(Fastq));
	buf=gzgets(fq,buf,1024*sizeof(char));
	if (!gzeof(fq)) {
		buf=gzgets(fq,buf,1024*sizeof(char));
		*(buf+strlen(buf)-1)='\0';	//remove the last \n in the buf
		line->seq = strdup(buf);
		buf=gzgets(fq,buf,1024*sizeof(char));
		buf=gzgets(fq,buf,1024*sizeof(char));
		*(buf+strlen(buf)-1)='\0';	//remove the last \n in the buf
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

double **load_fastq_file(const char *fq_file1,const char *fq_file2,
	unsigned long *total_reads_count,long long begin,HASHTBL **tbl,int **Quality,int **Length,int **Nucleotide) {
	char *buf = (char *)malloc(1024 * sizeof(char));
	double **PE_GC=(double **)calloc((strcmp(fq_file2,"")?2:1),sizeof(double *));
	int GCCountAlloc=ELECNT;
	Alloc2DArray(*Quality,128,MaxLen,*Nucleotide,5,MaxLen,*Length,MaxLen,*PE_GC,GCCountAlloc);
	HSIZE tblsiz=(HSIZE)ELECNT*1.34;
	*tbl = hashtbl_create(tblsiz, NULL);
	gzFile fastq1=creat_gzinfile(fq_file1);
	gzFile fastq2 = NULL;
	if (strcmp(fq_file2,"")) {
		fastq2=creat_gzinfile(fq_file2);
		Alloc2DArray(*(Quality+1),128,MaxLen,*(Nucleotide+1),5,MaxLen,*(Length+1),MaxLen,*(PE_GC+1),GCCountAlloc);
	}
	Fastq *line1=NULL,*line2=NULL;
	char *pair_seq=(char *)calloc(512 , sizeof(char));
	double total_GC=0;
	while (1){
		line1=readNextNode(fastq1,buf);
		if (line1==NULL) break;
		int Len1=0;
		STATSEQ(line1->seq,*(*PE_GC+(*total_reads_count)),*Nucleotide,Len1);
		total_GC+=*(*PE_GC+(*total_reads_count));
		Len1>75 ?memcpy(pair_seq,line1->seq,50):memcpy(pair_seq,line1->seq,Len1);
		AssignQuality(*Quality,line1->quality);
		(*(*Length+Len1-1))++;
		
		if (strcmp(fq_file2,"")) {
			line2=readNextNode(fastq2,buf);
			int Len2=0;
			STATSEQ(line2->seq,*(*(PE_GC+1)+(*total_reads_count)),*(Nucleotide+1),Len2);
			Len2>75 ? memcpy(pair_seq+50,line2->seq,50):memcpy(pair_seq+Len1,line2->seq,Len2);
			AssignQuality(*(Quality+1),line2->quality);
			(*(*(Length+1)+Len2-1))++;
		}
		void *data=hashtbl_get(*tbl,pair_seq);
		if (data==NULL){
			int *count=malloc(sizeof(int ));
			*count=1;
			hashtbl_insert(*tbl, pair_seq, (void *)count);
		}else{
			((*(int *)data))++;
		}
		free_Fastq(line1);
		if (strcmp(fq_file2,"")) free_Fastq(line2);
		(*total_reads_count)++;
		if (*total_reads_count>=GCCountAlloc){
			GCCountAlloc*=1.2;
			*PE_GC=realloc(*PE_GC,GCCountAlloc*sizeof(double));
			if (strcmp(fq_file2,"")) {
				*(PE_GC+1)=realloc(*(PE_GC+1),GCCountAlloc*sizeof(double));
			}
		}
		memset(pair_seq,0,512 * sizeof(char));
	}
	gzclose(fastq1);
	*PE_GC=realloc(*PE_GC,(*total_reads_count)*sizeof(double));
	if (strcmp(fq_file2,"")) {
		gzclose(fastq2);
		*(PE_GC+1)=realloc(*(PE_GC+1),(*total_reads_count)*sizeof(double));
	}
	free(pair_seq);
	free(buf);
	fprintf(stderr,"mean GC%% = %f%%\n",(double)total_GC/(*total_reads_count)*100);
	fprintf(stderr,"hash size: %ld\n",(unsigned long)(*tbl)->size);
	fprintf(stderr,"unique reads %d (%d/%ld= %.3f%% )\n",(int)(*tbl)->count,(int)(*tbl)->count,*total_reads_count,(double)(*tbl)->count/(*total_reads_count)*100);
	fprintf(stderr,"Finished load hash at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	return PE_GC;
}

int compare_hashed_data_count(const void *a, const void *b) {
	return *((int *)(*(ENTRY* const *)b)->data)-*((int *)(*(ENTRY* const *)a)->data);
}

void hash_data_delete(void *DATA) {
	ENTRY *node=(ENTRY *)DATA;
//	free(node->data);
	free(node);
}

void int2str(int c, int base, char *ret) {
	const char *tab = "0123456789abcdef";
	if (c == 0) ret[0] = '0', ret[1] = 0;
	else {
		int l=0, x=0, y=0;
		char buf[16];
		for (l = 0, x = c < 0? -c : c; x > 0; ) {
			buf[l++] = tab[x%base];
			x /= base;
//			/**********insert thousands**************/
//			if (x && (l-y)%3==0 ) {
//				y++;
//				buf[l++] = ',' ;
//			}
		}
		if (c < 0) buf[l++] = '-';
		for (x = l - 1, y = 0; x >= 0; --x) ret[y++] = buf[x];
		ret[y] = 0;
	}
}

SEXP qsort_hash_count(SEXP fq1,SEXP fq2){
	LOAD_HASHTBL(fq1,fq2);
	ENTRY **nodes=dump_hash_table(hashtbl);
	qsort(nodes, hashtbl->count, sizeof(ENTRY *), compare_hashed_data_count);
	unsigned long i, nprot=0;
	SEXP myint,mydouble,mydouble1,list,RQuality,RQuality2,RNucleotide,RNucleotide2,RLength,RLength2,LenName;
	PROTECT(myint = NEW_INTEGER(hashtbl->count));nprot++;
	int *p_myint = INTEGER_POINTER(myint);
	for (i=0;i<hashtbl->count ;i++ ) {
		p_myint[i] = *(int *)(nodes[i]->data);
	}
	PROTECT(LenName = allocVector(STRSXP,MaxLen));nprot++;
	char tmp[16]={0};
	for (i=0;i<MaxLen ;i++ ){
		int2str(i+1,10,tmp);
		SET_STRING_ELT(LenName,i,mkChar(tmp));
	}

	PROTECT(mydouble = NEW_NUMERIC(total_reads_num));nprot++;
	memmove(REAL(mydouble),*PE_GC,total_reads_num*sizeof(double));
	PROTECT(RQuality=allocMatrix(INTSXP,128,MaxLen));nprot++;
	memmove(INTEGER(RQuality),*Quality,128*MaxLen*sizeof(int));
	PROTECT(RNucleotide=allocMatrix(INTSXP,5,MaxLen));nprot++;
	memmove(INTEGER(RNucleotide),*Nucleotide,5*MaxLen*sizeof(int));
	PROTECT(RLength = allocVector(INTSXP,MaxLen));nprot++;
	memmove(INTEGER(RLength),*Length,MaxLen*sizeof(int));
	setAttrib(RLength, R_NamesSymbol, LenName);

	PROTECT(list = allocVector(VECSXP, strcmp(fq_file2,"")?9:5));nprot++;
	SET_VECTOR_ELT(list, 0, myint);
	SET_VECTOR_ELT(list, 1, mydouble);
	SET_VECTOR_ELT(list, 2, RQuality);
	SET_VECTOR_ELT(list, 3, RNucleotide);
	SET_VECTOR_ELT(list, 4, RLength);
	if (strcmp(fq_file2,"")) {
		PROTECT(mydouble1 = NEW_NUMERIC(total_reads_num));nprot++;
		memmove(REAL(mydouble1),*(PE_GC+1),total_reads_num*sizeof(double));
		PROTECT(RQuality2=allocMatrix(INTSXP,128,MaxLen));nprot++;
		memmove(INTEGER(RQuality2),*(Quality+1),128*MaxLen*sizeof(int));
		PROTECT(RNucleotide2=allocMatrix(INTSXP,5,MaxLen));nprot++;
		memmove(INTEGER(RNucleotide2),*(Nucleotide+1),5*MaxLen*sizeof(int));
		PROTECT(RLength2 = allocVector(INTSXP,MaxLen));nprot++;
		memmove(INTEGER(RLength2),*(Length+1),MaxLen*sizeof(int));
		setAttrib(RLength2, R_NamesSymbol, LenName);
		SET_VECTOR_ELT(list, 5, mydouble1);
		SET_VECTOR_ELT(list, 6, RQuality2);
		SET_VECTOR_ELT(list, 7, RNucleotide2);
		SET_VECTOR_ELT(list, 8, RLength2);
		
	}
	UNPROTECT(nprot);
	free(nodes);
	hashtbl_destroy(hashtbl,hash_data_delete);
	fprintf(stderr,"Finished at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	return list;
}

