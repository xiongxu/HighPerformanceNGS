#ifndef _IO_H
#define _IO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
//#include <err.h>
#include <fcntl.h>
#include <unistd.h>
/*
typedef pFILE FILE *;
#define OPEN_STREAM(TYPE,fileno,mode)	\
TYPE open_stream_##TYPE(char* filename,void(*f)(int , const char *)) {		\
	int fd ;																\
	if (strncmp(filename,"-", 1)==0) {										\
		fd = fileno;														\
	} else {																\
		fd = open(filename, O_CREAT | O_WRONLY | O_TRUNC, 0666 );			\
		if (fd==-1) err(1, "Failed to create output file (%s)", filename);	\
	}																		\
	TYPE out = f(fd,mode);													\
	return out;																\
}

#define CREATE_OUTFILE(TYPE)	\
TYPE creat_outfile_##TYPE(char *outfile,char *suffix) {
	char *out=(char *)calloc(strlen(outfile)+strlen(suffix)+2,sizeof(char));
	sprintf(out,"%s%s",outfile,suffix);
	if (!strncmp(outfile,"-",1)){
		TYPE fo=open_stream_##TYPE(out,fdopen);
	}
	else{
		TYPE fo=open_stream_##TYPE(out,gzdopen);
	}
	return fo;
}

OPEN_STREAM(pFILE,STDOUT_FILENO,"wb")
OPEN_STREAM(gzFile,STDOUT_FILENO,"wb")

CREATE_OUTFILE(gzFile)
CREATE_OUTFILE(pFILE)
*/
static inline FILE *fopen_input_stream(const char *filename);
static inline FILE *fopen_output_stream(const char *filename);
static inline FILE *fcreat_infile(const char *infile,const char *suffix);
static inline FILE *fcreat_outfile(const char *outfile,const char *suffix);
static inline gzFile open_output_stream(char* filename);
static inline gzFile creat_outfile(char *outfile,char *suffix);
static inline gzFile open_input_stream(const char *filename);

FILE *fopen_input_stream(const char *filename){
	int fd ;
	if (strncmp(filename,"-", 1)==0 || !strcmp(filename,"")) {
		fd = STDIN_FILENO;
	} else {
#if defined(__APPLE__) && defined(__MACH__) || defined(unix) || defined(linux)
		fd = open(filename, O_CREAT | O_RDONLY , 0666 );
#else
		fd = _open(filename, _O_CREAT|_O_BINARY | _O_RDONLY , 0666 );
#endif
		if (fd==-1) fprintf(stderr, "Failed to create input file (%s)", filename);
	}
	FILE *in = fdopen(fd,"rb");
	return in;
}

FILE *fopen_output_stream(const char *filename) {
	int fd ;
	if (strncmp(filename,"-", 1)==0 || !strcmp(filename,"")) {
		fd = STDOUT_FILENO;
	} else {
#if defined(__APPLE__) && defined(__MACH__) || defined(unix) || defined(linux)
		fd = open(filename, O_CREAT | O_WRONLY | O_TRUNC, 0666 );
#else
		fd = _open(filename, _O_CREAT|_O_BINARY | _O_WRONLY | _O_APPEND |_O_TRUNC, 0666 );
#endif
		if (fd==-1) fprintf(stderr, "Failed to create output file (%s)", filename);
	}
	FILE *out = fdopen(fd,"wb");
	return out;
}

FILE *fcreat_infile(const char *infile,const char *suffix) {
	char in[strlen(infile)+strlen(suffix)+2];
	sprintf(in,"%s%s",infile,suffix);
	FILE *fi=fopen_output_stream(in);
	return fi;
}

FILE *fcreat_outfile(const char *outfile,const char *suffix) {
	char *out=(char *)calloc(strlen(outfile)+strlen(suffix)+2,sizeof(char));
	sprintf(out,"%s%s",outfile,suffix);
	FILE *fo=fopen_output_stream(out);
	return fo;
}

gzFile open_output_stream(char* filename) {
	int fd ;
	if (strncmp(filename,"-", 1)==0 || !strcmp(filename,"")) {
		fd = STDOUT_FILENO;
	} else {
#if defined(__APPLE__) && defined(__MACH__) || defined(unix) || defined(linux)
		fd = open(filename, O_CREAT | O_WRONLY | O_TRUNC, 0666 );
#else
		fd = _open(filename, _O_CREAT|_O_BINARY | _O_WRONLY | _O_APPEND |_O_TRUNC, 0666 );
#endif
		if (fd==-1) fprintf(stderr, "Failed to create output file (%s)", filename);
	}
	gzFile out = gzdopen(fd,"wb");
	return out;
}

gzFile creat_outfile(char *outfile,char *suffix) {
	char *out=(char *)calloc(strlen(outfile)+strlen(suffix)+2,sizeof(char));
	sprintf(out,"%s%s",outfile,suffix);
	gzFile fo=open_output_stream(out);
	return fo;
}

gzFile open_input_stream(const char *filename){
	int fd ;
	if (strncmp(filename,"-", 1)==0 || !strcmp(filename,"")) {
		fd = STDIN_FILENO;
	} else {
#if defined(__APPLE__) && defined(__MACH__) || defined(unix) || defined(linux)
		fd = open(filename, O_CREAT | O_RDONLY, 0666 );
#else
		fd = _open(filename, _O_CREAT|_O_BINARY | _O_RDONLY , 0666 );
#endif
		if (fd==-1) fprintf(stderr, "Failed to create input file (%s)", filename);
	}
	gzFile in = gzdopen(fd,"rb");
	return in;
}
#endif
