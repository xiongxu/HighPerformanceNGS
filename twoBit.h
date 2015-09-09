#ifndef _TWOBIT_H
#define _TWOBIT_H
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include "hiredis.h"
#include "sds.h"

#define UBYTE unsigned char   /* Wants to be unsigned 8 bits. */
#define BYTE signed char      /* Wants to be signed 8 bits. */
#define UWORD unsigned short  /* Wants to be unsigned 16 bits. */
#define WORD short        /* Wants to be signed 16 bits. */
#define bits64 unsigned long long  /* Wants to be unsigned 64 bits. */
#define bits32 unsigned       /* Wants to be unsigned 32 bits. */
#define bits16 unsigned short /* Wants to be unsigned 16 bits. */
#define bits8 unsigned char   /* Wants to be unsigned 8 bits. */
#define signed32 int          /* Wants to be signed 32 bits. */

#define T_BASE_VAL 0
#define U_BASE_VAL 0
#define C_BASE_VAL 1
#define A_BASE_VAL 2
#define G_BASE_VAL 3
#define N_BASE_VAL 4   /* Used in 1/2 byte representation. */

#define packed4Size(unpackedSize) ((unpackedSize + 3) >> 2)
#define packed3Size(unpackedSize) ((unpackedSize + 2) / 3)

static inline int seq2twoBit(char *seq,UBYTE **packed_seq);
static inline void unpackDna4(UBYTE *tiles, int byteCount, char *out);
static inline UBYTE packDna4(char *in);
static inline char *twoBit2Seq(UBYTE *packed,int fragStart, int fragEnd);
static inline void initNtVal(void);

static inline UBYTE packDna3(char *in);
static inline int seq2packDna3(char *seq,UBYTE **packed_seq);
static inline char *twoBit2Seq3(UBYTE *packed,int fragStart, int fragEnd);
static inline void int2str(int c, int base, char *ret);
static inline UBYTE *int2char(int pos,UBYTE *value);
static inline int char2int(UBYTE *string,int *pos);
static inline uint32_t sds2int(const sds string);
static inline sds int2sds(uint32_t pos);

static inline char *sds2seq(sds packed,int fragStart, int fragEnd);
static inline char *unpackSds(sds packedSeq);
static inline sds packSds(char *seq);
static inline sds seq2sds(char *seq);

extern int ntValNoN[256];
extern char valToNt[5];
extern char revNt[256];

void initNtVal(void){
	int i;
//	fprintf(stderr,"sizeof(UBYTE): %d\n",(int)sizeof(UBYTE));
	for (i=0; i<256; i++){
		ntValNoN[i] = T_BASE_VAL;
		revNt[i]=0;
	}
	ntValNoN['t'] = ntValNoN['T'] = T_BASE_VAL;
	ntValNoN['u'] = ntValNoN['U'] = U_BASE_VAL;
	ntValNoN['c'] = ntValNoN['C'] = C_BASE_VAL;
	ntValNoN['a'] = ntValNoN['A'] = A_BASE_VAL;
	ntValNoN['g'] = ntValNoN['G'] = G_BASE_VAL;
	valToNt[T_BASE_VAL]= revNt['A']=revNt['a'] = 'T';
	valToNt[C_BASE_VAL]= revNt['G']=revNt['g'] = 'C';
	valToNt[A_BASE_VAL]= revNt['T']=revNt['t'] = 'A';
	valToNt[G_BASE_VAL]= revNt['C']=revNt['c'] = 'G';
	valToNt[N_BASE_VAL]= revNt['N']=revNt['.'] = 'N';
}

char *twoBit2Seq(UBYTE *packed,int fragStart, int fragEnd){
	int i,j;
	int packByteCount, packedStart, packedEnd, remainder, midStart, midEnd;
	char *dna=(char *)calloc(fragEnd - fragStart+1,sizeof(char));
	char *pt_dna=dna;

	/* Skip to bits we need and read them in. */
	packedStart = (fragStart>>2);
	packedEnd = ((fragEnd+3)>>2);
	packByteCount = packedEnd - packedStart;
	packed+=packedStart;

	/* Handle case where everything is in one packed byte */
	if (packByteCount == 1){
		int pOff = (packedStart<<2);
		int pStart = fragStart - pOff;
		int pEnd = fragEnd - pOff;
		UBYTE partial = *packed;
		assert(pEnd <= 4);
		assert(pStart >= 0);
		for (i=pStart; i<pEnd; ++i)
			*pt_dna++ = valToNt[(partial >> (6-i-i)) & 3];
	}
	else{
		/* Handle partial first packed byte. */
		midStart = fragStart;
		remainder = (fragStart&3);
		if (remainder > 0){
			UBYTE partial = *packed++;
			int partCount = 4 - remainder;
			for (i=partCount-1; i>=0; --i){
				pt_dna[i] = valToNt[partial&3];
				partial >>= 2;
			}
			midStart += partCount;
			pt_dna += partCount;
		}

		/* Handle middle bytes. */
		remainder = fragEnd&3;
		midEnd = fragEnd - remainder;
		for (i=midStart; i<midEnd; i += 4){
			UBYTE b = *packed++;
			for (j=3; j>=0; --j){
				pt_dna[j] = valToNt[b & 0x3];
				b >>= 2;
			}
			pt_dna += 4;
		}

		if (remainder >0){
			UBYTE part = *packed;
			part >>= (8-remainder-remainder);
			for (i=remainder-1; i>=0; --i){
				pt_dna[i] = valToNt[part&3];
				part >>= 2;
			}
		}
	}
	return dna;
}

UBYTE packDna4(char *in){
	UBYTE out = 0;
	short i;
	for (i=0;i<4 ;i++){
		out = out <<2 | ntValNoN[(short)*(in+i)];
	}
	return out;
}

sds packSds(char *seq){
	int seq_size=strlen(seq);
	int codeLen = (int)packed4Size(seq_size);
	sds packed_seq=sdsnewlen(NULL,codeLen);
	int i=0,j=0;
	for (i=0; i<seq_size; i += 4){
		packed_seq[j++] = (char)packDna4(seq+i);
	}
	return packed_seq;
}

char *unpackSds(sds packedSeq){
	int i, j;
	char *out=(char *)calloc(37,sizeof(char));
	int unpackSdsLen=sdslen(packedSeq);
	for (i=0; i<unpackSdsLen; ++i){
		for (j=0; j<4; j++){
			out[i*4+j] = valToNt[packedSeq[i] >> (3-j)*2 & 0x3];
		}
	}
	return out;
}

sds seq2sds(char *seq){
	size_t seq_size=strlen(seq);
	int codeLen = (int)packed4Size(seq_size);
	sds packed_seq=sdsnewlen(NULL,codeLen);
	char last4[4];
	int i=0,j=0,end = seq_size - 4;
	for (i=0; i<end; i += 4){
		packed_seq[j++] = (char)packDna4(seq+i);
	}
	/* handle the last 1/2/3/4 base */
	memset(last4,'T',4);
	memcpy(last4, seq+i, seq_size-i);
	packed_seq[j++] = (char)packDna4(last4);
	return packed_seq;
}

char *sds2seq(sds packed,int fragStart, int fragEnd){
	int i,j;
	int packByteCount, packedStart, packedEnd, remainder, midStart, midEnd;
	char *dna=(char *)calloc(fragEnd - fragStart+1,sizeof(char));
	char *pt_dna=dna;

	/* Skip to bits we need and read them in. */
	packedStart = (fragStart>>2);
	packedEnd = ((fragEnd+3)>>2);
	packByteCount = packedEnd - packedStart;
	packed+=packedStart;

	/* Handle case where everything is in one packed byte */
	if (packByteCount == 1){
		int pOff = (packedStart<<2);
		int pStart = fragStart - pOff;
		int pEnd = fragEnd - pOff;
		UBYTE partial = (UBYTE)*packed;
		assert(pEnd <= 4);
		assert(pStart >= 0);
		for (i=pStart; i<pEnd; ++i)
			*pt_dna++ = valToNt[(partial >> (6-i-i)) & 3];
	}
	else{
		/* Handle partial first packed byte. */
		midStart = fragStart;
		remainder = (fragStart&3);
		if (remainder > 0){
			UBYTE partial = (UBYTE)*packed++;
			int partCount = 4 - remainder;
			for (i=partCount-1; i>=0; --i){
				pt_dna[i] = valToNt[partial&3];
				partial >>= 2;
			}
			midStart += partCount;
			pt_dna += partCount;
		}

		/* Handle middle bytes. */
		remainder = fragEnd&3;
		midEnd = fragEnd - remainder;
		for (i=midStart; i<midEnd; i += 4){
			UBYTE b = (UBYTE)*packed++;
			for (j=3; j>=0; --j){
				pt_dna[j] = valToNt[b & 0x3];
				b >>= 2;
			}
			pt_dna += 4;
		}

		if (remainder >0){
			UBYTE part = (UBYTE)*packed;
			part >>= (8-remainder-remainder);
			for (i=remainder-1; i>=0; --i){
				pt_dna[i] = valToNt[part&3];
				part >>= 2;
			}
		}
	}
	return dna;
}

UBYTE packDna3(char *in){
	UBYTE out = 0;
	int count = 3;
	int bVal;
	while (--count >= 0) {
		bVal = ntValNoN[(int)*in++];
		out <<= 2;
		out += bVal;
	}
	return out+33;
}

void unpackDna4(UBYTE *tiles, int byteCount, char *out){
	int i, j;
	UBYTE tile;

	for (i=0; i<byteCount; ++i){
		tile = tiles[i];
		for (j=3; j>=0; --j){
			out[j] = valToNt[tile & 0x3];
			tile >>= 2;
		}
		out += 4;
	}
}


int seq2twoBit(char *seq,UBYTE **packed_seq){
	size_t seq_size=strlen(seq);
	int ubyteSize = (int)packed4Size(seq_size);
	*packed_seq=(UBYTE *)calloc(ubyteSize,sizeof(UBYTE));
	char last4[4];
	int i=0,j=0,end = seq_size - 4;
	for (i=0; i<end; i += 4){
		*(*packed_seq+j++) = packDna4(seq+i);
	}
	/* handle the last 1/2/3/4 base */
	memset(last4,'T',4);
	memcpy(last4, seq+i, seq_size-i);
	*(*packed_seq+j++) = packDna4(last4);
	return j;
}

int seq2packDna3(char *seq,UBYTE **packed_seq){
	size_t seq_size=strlen(seq);
	int ubyteSize = (int)packed3Size(seq_size);
	*packed_seq=(UBYTE *)calloc(ubyteSize,sizeof(UBYTE));
	char last3[3];
	int i=0,j=0,end = seq_size - 3;
	for (i=0,j=0; i<end; i += 3,j++){
		*(*packed_seq+j) = packDna3(seq+i);
	}
	/* handle the last 1/2/3/4 base */
	memset(last3,'T',3);
	memcpy(last3, seq+i, seq_size-i);
	*(*packed_seq+j++) = packDna3(last3);
	return j;
}


char *twoBit2Seq3(UBYTE *packed,int fragStart, int fragEnd){
	int i,j;
	int packByteCount, packedStart, packedEnd, remainder, midStart, midEnd;
	char *dna=(char *)calloc(fragEnd - fragStart+1,sizeof(char));
	char *pt_dna=dna;

	/* Skip to bits we need and read them in. */
	packedStart = (fragStart/3);
	packedEnd = ((fragEnd+2)/3);
	packByteCount = packedEnd - packedStart;
	packed+=packedStart;

	/* Handle case where everything is in one packed byte */
	if (packByteCount == 1){
		int pOff = (packedStart/3);
		int pStart = fragStart - pOff;
		int pEnd = fragEnd - pOff;
		UBYTE partial = *packed-33;
		assert(pEnd <= 3);
		assert(pStart >= 0);
		for (i=pStart; i<pEnd; ++i)
			*pt_dna++ = valToNt[(partial >> (4-i-i)) & 3];
	}
	else{
		/* Handle partial first packed byte. */
		midStart = fragStart;
		remainder = (fragStart&3);
		if (remainder > 0){
			UBYTE partial = *packed++ -33;
			int partCount = 3 - remainder;
			for (i=partCount-1; i>=0; --i){
				pt_dna[i] = valToNt[partial&3];
				partial >>= 2;
			}
			midStart += partCount;
			pt_dna += partCount;
		}

		/* Handle middle bytes. */
		remainder = fragEnd%3;
		midEnd = fragEnd - remainder;
		for (i=midStart; i<midEnd; i += 3){
			UBYTE b = *packed++-33;
			for (j=2; j>=0; --j){
				pt_dna[j] = valToNt[b&0x3];
				b >>= 2;
			}
			pt_dna += 3;
		}

		if (remainder >0){
			UBYTE part = *packed-33;
			part >>= (6-remainder-remainder);
			for (i=remainder-1; i>=0; --i){
				pt_dna[i] = valToNt[part&3];
				part >>= 2;
			}
		}
	}
	return dna;
}

/*
int str2int(char *str){
	int temp = 0;
	const char *ptr = str;
	if (*str == '-' || *str == '+') str++;
	while(*str != 0){
		if ((*str < '0') || (*str > '9')){
			break;
		}
		temp = temp * 10 + (*str - '0');
		str++;
	}
	if (*ptr == '-') temp = -temp;
	return temp;
}
*/

void int2str(int c, int base, char *ret)
{
	const char *tab = "0123456789abcdef";
	if (c == 0) ret[0] = '0', ret[1] = 0;
	else {
		int l=0, x=0, y=0;
		char buf[16];
		for (l = 0, x = c < 0? -c : c; x > 0; ) {
			buf[l++] = tab[x%base];
			x /= base;
			/**********insert thousands**************/
			if (x && (l-y)%3==0 ) {
				y++;
				buf[l++] = ',' ;
			}
		}
		if (c < 0) buf[l++] = '-';
		for (x = l - 1, y = 0; x >= 0; --x) ret[y++] = buf[x];
		ret[y] = 0;
	}
}

UBYTE *int2char(int pos,UBYTE *value){
	int i;
	for (i=3;i>=0 ;i-- )
		value[i]= (pos>>(i*7) & 0x7F) +1;
	value[4]=0;
	return value;
}

int char2int(UBYTE *string,int *pos){
	int i;
	*pos=0;
	for (i=3;i>=0 ;i-- ){
		*pos<<=7;
		*pos+=(string[i]-1);
	}
	return *pos;
}

sds int2sds(uint32_t pos){
	short codeLen=1;
	if (pos>16777215){
		codeLen=4;
	}else if(pos>65535){
		codeLen=3;
	}else if(pos>255){
		codeLen=2;
	}
	sds value=sdsnewlen(NULL,codeLen);
	short i;
	for (i=codeLen-1;i>=0 ;i-- )
		value[i]= pos>>(i*8) & 0xFF;
	return value;
}

uint32_t sds2int(const sds string){
	uint32_t pos=0;
	short i;
	for (i=sdslen(string)-1; i>=0 ;i-- ){
		pos=pos << 8 | (unsigned char)string[i];
	}
	return pos;
}


#ifdef TWOBIT_TEST
static inline void twobit_test(void);
void twobit_test(void) {
	int i,pos=0;
	for (i=0;i<100 ;i++ ){
		uint32_t random=1+(uint32_t)(100000.0*rand()/(RAND_MAX+1.0));
		fprintf(stdout,"%d\t",random);
		sds intchar=int2sds(random);
		fprintf(stdout,"%s\t",sdscatrepr(sdsempty(),intchar,sdslen(intchar)));
		pos=sds2int(intchar);
		fprintf(stdout,"%d\n",(int)pos);
		sdsfree(intchar);
	}
	char *seq="AAGTATCAAGTGAGTAATATGATGGGAAGACTTTTA";
	sds packedSeq=packSds(seq);
	fprintf(stdout,"%s\n",sdscatrepr(sdsempty(),packedSeq,sdslen(packedSeq)));
	char *unpackedSeq=unpackSds(packedSeq);
	fprintf(stdout,"%s\n%s\n%d\n",seq,unpackedSeq,strcmp(seq,unpackedSeq));
}
#endif
#endif
