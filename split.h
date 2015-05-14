#ifndef _SSPLIT_H
#define _SSPLIT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#define printSplitString(string,string_index)	\
do{													\
	int index;										\
	for (index=0;index<string_index ;index++ ){		\
		fprintf(stdout,"%s\n",*(string+index));		\
	}												\
}while (0)

#define freeSplitString(string,string_index)	\
do{												\
	int index;									\
	for (index=0;index<string_index ;index++ ){	\
		free(*(string+index));					\
	}											\
	free(string);								\
}while (0)

static inline int str2int(const char *str);
static inline int split(char **string,const char *str,const char sep);

int str2int(const char *str){
	int temp = 0;
	const char *ptr = str;
	if (*str == '-' || *str == '+') str++;
	while(*str != 0){
//		fprintf(stderr,"%c\n",*str);
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

int split(char **string,const char *str,const char sep){
	int i=0,t=0,string_index=0;
	while(*(str+t)) {
		if (*(str+t)==sep ) {
			if (i){
				*(string+string_index)=(char *)calloc(i+1,sizeof(char));
				memcpy(*(string+string_index),str+(t-i),i);
				*(*(string+string_index++)+i)=0;
			}
			i=0;
		}else{
			 i++;
		}
		t++;
	}
	if (i) {
		*(string+string_index)=(char *)calloc(i+1,sizeof(char));
		memcpy(*(string+string_index),str+(t-i),i);
		*(*(string+string_index++)+i)=0;
	}
	string=(char **)realloc(string,string_index*sizeof(char *));
	return string_index;
}

double str2double(const char *str){
	if (!strcmp(str,"NA")) {
		return 0;
	}else{
		return strtod(str, NULL);
	}
}

int split2double(double *string,char *str,char **chr,char sep,char temp[32]){
	unsigned int i=0,t=0,string_index=0,flag=0;
//	char *temp=(char *)calloc(32,sizeof(char));
	while(*(str+t) != 0) {
		if (*(str+t)==sep ) {
			if (i){ 
				strncpy(temp,str+(t-i),i);
				*(temp+i)=0;
				if (flag) {
					string[string_index++]=str2double(temp);
				}else{
					*chr=strdup(temp);
				}
			}
			flag++;
			i=0;
		}else{
			i++;
		}
		t++;
	}
	if (i) {
		strncpy(temp,str+(t-i),i);
		*(temp+i)=0;
		if (flag) {
			string[string_index++]=str2double(temp);
		}else{
			*chr=strdup(temp);
		}
	}
//	free(temp);
	return string_index;
}

#ifdef DSPLIT
static inline void test_split(int alloc);
void test_split(int alloc) {
	int string_index;
	const char* seq="AA CCC  GGGG TTTTTT   N";
	char **splitSeq=(char **)malloc(alloc*sizeof(char *));
	string_index=split(splitSeq,seq,' ',alloc);
	fprintf(stderr,"%d\n",string_index);
	printSplitString(splitSeq,string_index);
	freeSplitString(splitSeq,string_index);
}
#endif

#endif