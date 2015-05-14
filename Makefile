CC=		gcc
CFLAGS=		-O3 -Wall -W -Wstrict-prototypes -Wwrite-strings -g -ggdb
LDFLAGS=	-lz 
DFLAGS=		
INCLUDES=	-I$(CURDIR)/zlib-1.2.8/include
LIBPATH=	-L$(CURDIR)/zlib-1.2.8/lib
LIBNAME=	liblist
STLIBSUFFIX=a
STLIBNAME=$(LIBNAME).$(STLIBSUFFIX)
STLIB_MAKE_CMD=ar rcs $(STLIBNAME)

.PHONY: all clean
all:zlib-1.2.8 hiredisDir gzfastq_sort gzfastq_sort_list pick_pair fastq_count fastq_count_kthread fastq_trim gzfastq_uniq gzfastq_sample bam2depth bam_sliding_count $(STLIBNAME) libhashtbl.a
SUBDIRS = `find -maxdepth 1 -type d | sed "1d"|grep -P "hiredis|samtools-0.1.19"`
zlib-1.2.8:zlib-1.2.8.tar.gz
	tar -zxvf zlib-1.2.8.tar.gz && mv zlib-1.2.8 zlib;\
	cd zlib;\
	test -d $(CURDIR)/zlib-1.2.8 || mkdir -p $(CURDIR)/zlib-1.2.8;\
	./configure --prefix=$(CURDIR)/zlib-1.2.8;\
	make;\
	make install;\
	cd $(CURDIR) && rm -rf zlib
	
hiredisDir:
	wdir=`pwd`; \
	for subdir in $(SUBDIRS); do \
		cd  $$subdir;\
		$(MAKE)  ;\
		cd $$wdir;\
	done
#.c.o:
#	 $(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
list.o:list.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $< -o $@
hashtbl.o:hashtbl.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $< -o $@

pick_pair.o:pick_pair.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
pick_pair:pick_pair.o
	$(CC) $(CFLAGS) ${DFLAGS} $< -o $@ $(LIBPATH) $(LDFLAGS)
gzfastq_sort.o:gzfastq_sort.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
gzfastq_sort:gzfastq_sort.o
	$(CC) $(CFLAGS) ${DFLAGS} $< -o $@ $(LIBPATH) $(LDFLAGS)
klib/kthread.o:klib/kthread.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $< -o $@
fastq_count_kthread.o:fastq_count_kthread.c 
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@ 
fastq_count_kthread:fastq_count_kthread.o klib/kthread.o
	$(CC) $(CFLAGS) $(DFLAGS) $^ -o $@ $(LIBPATH) $(LDFLAGS) -lpthread
fastq_count.o:fastq_count.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
fastq_count:fastq_count.o
	$(CC) $(CFLAGS) ${DFLAGS} $< -o $@ $(LIBPATH) $(LDFLAGS) -lpthread

fastq_trim.o:fastq_trim.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
fastq_trim:fastq_trim.o
	$(CC) $(CFLAGS) ${DFLAGS} $< -o $@ $(LIBPATH) $(LDFLAGS)
gzfastq_uniq.o:gzfastq_uniq.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -I./hiredis $< -o $@
gzfastq_uniq:gzfastq_uniq.o
	$(CC) $(CFLAGS) ${DFLAGS} $< -o $@ $(LIBPATH) -L./hiredis $(LDFLAGS) -lhiredis

gzfastq_sample.o:gzfastq_sample.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -I./klib -I./fastq-tools-0.7/src $< -o $@
gzfastq_sample:gzfastq_sample.o
	$(CC) $(CFLAGS) ${DFLAGS} $< -o $@ $(LIBPATH) -L./fastq-tools-0.7/src $(LDFLAGS) -lrng

gzfastq_sort_list.o:gzfastq_sort_list.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -I. $< -o $@
gzfastq_sort_list:gzfastq_sort_list.o list.o
	$(CC) $(CFLAGS) ${DFLAGS} $^ -o $@ $(LIBPATH) $(LDFLAGS)

bam_sliding_count.o:bam_sliding_count.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -I./samtools-0.1.19 $< -o $@
bam_sliding_count:bam_sliding_count.o
	$(CC) $(CFLAGS) ${DFLAGS} $< -o $@ $(LIBPATH) -L./samtools-0.1.19 -lbam -lpthread -lgd -lpng -lfreetype -lm $(LDFLAGS)

bam2depth.o:bam2depth.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -I./samtools-0.1.19 -I. $< -o $@
bam2depth:bam2depth.o hashtbl.o
	$(CC) $(CFLAGS) $(DFLAGS) $^ -o $@ $(LIBPATH) -L./samtools-0.1.19 -L. -lbam -lpthread -lm  $(LDFLAGS) 

$(STLIBNAME):list.o
	$(STLIB_MAKE_CMD) $<

libhashtbl.a:hashtbl.o
	 $(AR) -csru $@ $<

clean:
	rm -rf hiredisDir gzfastq_sort gzfastq_sort_list pick_pair fastq_count fastq_count_kthread fastq_trim gzfastq_uniq gzfastq_sample bam_sliding_count bam2depth $(STLIBNAME) libhashtbl.a *.o; \
	wdir=`pwd`; \
	cd ./hiredis;\
	$(MAKE) clean;\
	cd $$wdir;
