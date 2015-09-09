CC=		gcc
CFLAGS=		-O3 -Wall -W -Wstrict-prototypes -Wwrite-strings -g -ggdb
LDFLAGS=	-lz 
DFLAGS=		
INCLUDES=	-I. -I$(CURDIR)/zlib-1.2.8/include -I/usr/local/include
LIBPATH=	-L. -L$(CURDIR)/zlib-1.2.8/lib -L/usr/local/lib
LIBNAME=	liblist
STLIBSUFFIX=a
STLIBNAME=$(LIBNAME).$(STLIBSUFFIX)
STLIB_MAKE_CMD=ar rcs $(STLIBNAME)

.PHONY: all clean
all:  hiredisDir gzfastq_sort gzfastq_sort_list pick_pair fastq_count fastq2twobit twoBit2seq gzfastq_uniq_sort fastq_count_kthread fastq_trim gzfastq_uniq gzfastq_sample bam2depth bam2wig bam_sliding_count $(STLIBNAME) Rgzfastq_uniq.so
SUBDIRS = `find $(CURDIR) -maxdepth 1 -type d | sed "1d"|grep -E "hiredis|samtools-0.1.19"`
$(CURDIR)/zlib-1.2.8:zlib-1.2.8.tar.gz
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
		$(MAKE) INCLUDES="$(INCLUDES)" LIBPATH="$(LIBPATH)" ;\
		cd $$wdir;\
	done

#.c.o:
#	 $(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
list.o:list.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $< -o $@
hashtbl.o:hashtbl.c
	$(CC) -fPIC -c $(CFLAGS) $(DFLAGS) $< -o $@

pick_pair.o:pick_pair.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
pick_pair:pick_pair.o $(CURDIR)/zlib-1.2.8
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
gzfastq_uniq_sort.o: gzfastq_uniq_sort.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
gzfastq_uniq_sort : gzfastq_uniq_sort.o libhashtbl.a
	$(CC) $(CFLAGS) ${DFLAGS} $< -o $@ $(LIBPATH) $(LDFLAGS) -lhashtbl
gzfastq_uniq.o:gzfastq_uniq.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -I./hiredis $< -o $@
gzfastq_uniq:gzfastq_uniq.o
	$(CC) $(CFLAGS) ${DFLAGS} $< -o $@ $(LIBPATH) -L./hiredis $(LDFLAGS) -lhiredis
fastq2twobit.o:fastq2twobit.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -I./hiredis $< -o $@
fastq2twobit:fastq2twobit.o list.o
	$(CC) $(CFLAGS) ${DFLAGS} $^ -o $@ $(LIBPATH) -L./hiredis $(LDFLAGS) -lhiredis
twoBit2seq.o:twoBit2seq.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -I./hiredis $< -o $@
twoBit2seq:twoBit2seq.o
	$(CC) $(CFLAGS) ${DFLAGS} $^ -o $@ $(LIBPATH) -L./hiredis $(LDFLAGS) -lhiredis

$(CURDIR)/fastq-tools-0.7:fastq-tools-0.7.tar.gz
	tar -zxvf fastq-tools-0.7.tar.gz ;\
	cd fastq-tools-0.7 ;\
	./configure ;\
	$(MAKE) ;\
	cd src && ar -csru librng.a common.o rng.o ;\
	cd $(CURDIR)

gzfastq_sample.o : gzfastq_sample.c $(CURDIR)/fastq-tools-0.7
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -I./klib -I./fastq-tools-0.7/src $< -o $@
gzfastq_sample : gzfastq_sample.o $(CURDIR)/fastq-tools-0.7/src/librng.a
	$(CC) $(CFLAGS) ${DFLAGS} $< -o $@ $(LIBPATH) -L./fastq-tools-0.7/src $(LDFLAGS) -lrng

gzfastq_sort_list.o:gzfastq_sort_list.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
gzfastq_sort_list:gzfastq_sort_list.o list.o
	$(CC) $(CFLAGS) ${DFLAGS} $^ -o $@ $(LIBPATH) $(LDFLAGS)

$(CURDIR)/samtools-0.1.19: $(CURDIR)/samtools-0.1.19.tar.bz2 $(CURDIR)/zlib-1.2.8
	tar -jxvf samtools-0.1.19.tar.bz2;\
	cd $(CURDIR)/samtools-0.1.19;\
	make INCLUDES="$(INCLUDES)" LIBPATH="$(LIBPATH)";\
	rm *.o misc/*.o bcftools/*.o;\
	cd $(CURDIR)

bam_sliding_count.o : bam_sliding_count.c $(CURDIR)/samtools-0.1.19
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -I./samtools-0.1.19 $< -o $@
bam_sliding_count:bam_sliding_count.o $(CURDIR)/samtools-0.1.19/libbam.a
	$(CC) $(CFLAGS) ${DFLAGS} $< -o $@ $(LIBPATH) -L./samtools-0.1.19 -lbam -lpthread -lgd -lpng -lfreetype -lm $(LDFLAGS)

bam2depth.o : bam2depth.c $(CURDIR)/samtools-0.1.19
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -I./samtools-0.1.19 $< -o $@
bam2depth : bam2depth.o hashtbl.o $(CURDIR)/samtools-0.1.19/libbam.a
	$(CC) $(CFLAGS) $(DFLAGS) $^ -o $@ $(LIBPATH) -L./samtools-0.1.19 -lbam -lpthread -lm  $(LDFLAGS) 

bam2wig.o:bam2wig.c $(CURDIR)/samtools-0.1.19
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -I./samtools-0.1.19 -I. $< -o $@
bam2wig:bam2wig.o hashtbl.o $(CURDIR)/samtools-0.1.19/libbam.a
	$(CC) $(CFLAGS) $(DFLAGS) $^ -o $@ $(LIBPATH) -L./samtools-0.1.19 -lbam -lpthread -lm  $(LDFLAGS)

$(STLIBNAME):list.o
	$(STLIB_MAKE_CMD) $<

libhashtbl.a:hashtbl.o
	 $(AR) -csru $@ $<

Rgzfastq_uniq.so:Rgzfastq_uniq.c libhashtbl.a
	gcc -shared -std=gnu99 -fPIC -g -O2 $< -o $@ $(INCLUDES) -I/usr/share/R/include -DNDEBUG  -L/usr/local/lib64 -L/usr/lib/R/lib -L. -lR -lhashtbl $(LDFLAGS)
clean:
	rm -rf fastq-tools-0.7 hiredisDir gzfastq_sort gzfastq_sort_list pick_pair fastq_count fastq_count_kthread fastq_trim fastq2twobit twoBit2seq gzfastq_uniq gzfastq_sample bam_sliding_count bam2depth bam2wig gzfastq_uniq_sort $(STLIBNAME) samtools-0.1.19 zlib-1.2.8  libhashtbl.a Rgzfastq_uniq.so *.o fastq-tools-0.7/src/librng.a ; \
	wdir=`pwd`; \
	cd ./hiredis;\
	$(MAKE) clean;\
	cd $$wdir;
