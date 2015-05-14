#include <stdio.h>
#include "sam.h"

typedef struct {
  uint32_t ltid;
  int      lstart,lcoverage,lpos,beg,end;
  samfile_t *in;
} tmpstruct_t;

// callback for bam_fetch()
static int fetch_func(const bam1_t *b, void *data)
{
	bam_plbuf_t *buf = (bam_plbuf_t*)data;
	bam_plbuf_push(b, buf);
	return 0;
}
// callback for bam_plbuf_init()
static int pileup_func(uint32_t tid, uint32_t position, int n, const bam_pileup1_t *pl, void *data)
{
  tmpstruct_t *tmp = (tmpstruct_t*)data;
  int pos          = (int)position;
  int coverage     = n;
  int i;
  for (i=0;i<n;i++)
    if (pl[i].is_del) coverage--;

  if (tmp->ltid != tid || tmp->lcoverage != coverage || pos > tmp->lpos+1) {
    if (tmp->lpos > 0 && tmp->lcoverage > 0)
      printf("%s\t%d\t%d\t%d\n", tmp->in->header->target_name[tmp->ltid], tmp->lstart,tmp->lpos+1, tmp->lcoverage);
    tmp->ltid       = tid;
    tmp->lstart     = pos;
    tmp->lcoverage  = coverage;
  }
  tmp->lpos = pos;
  return 0;
}

int main(int argc, char *argv[])
{
	tmpstruct_t tmp;
	if (argc == 1) {
		fprintf(stderr, "Usage: bam2bedgraph <in.bam> [region]\n\nCreate a BEDGraph file of genomic coverage. BAM file must be sorted.\n");
		return 1;
	}
	tmp.beg = 0; tmp.end = 0x7fffffff;
	tmp.lstart    = 0; 
	tmp.lcoverage = 0;
	tmp.ltid      = 0;
	tmp.lpos      = 0;

	tmp.in = samopen(argv[1], "rb", 0);
	if (tmp.in == 0) {
		fprintf(stderr, "Fail to open BAM file %s\n", argv[1]);
		return 1;
	}
	if (argc == 2) { // if a region is not specified
		sampileup(tmp.in, -1, pileup_func, &tmp);
	} else {
		int ref;
		bam_index_t *idx;
		bam_plbuf_t *buf;
		idx = bam_index_load(argv[1]); // load BAM index
		if (idx == 0) {
			fprintf(stderr, "BAM indexing file is not available.\n");
			return 1;
		}
		bam_parse_region(tmp.in->header, argv[2], &ref, &tmp.beg, &tmp.end); // parse the region
		if (ref < 0) {
			fprintf(stderr, "Invalid region %s\n", argv[2]);
			return 1;
		}
		buf = bam_plbuf_init(pileup_func, &tmp); // initialize pileup
		bam_fetch(tmp.in->x.bam, idx, ref, tmp.beg, tmp.end, buf, fetch_func);
		bam_plbuf_push(0, buf); // finalize pileup
		bam_index_destroy(idx);
		bam_plbuf_destroy(buf);
	}
	printf("%s\t%d\t%d\t%d\n", tmp.in->header->target_name[tmp.ltid], tmp.lstart,tmp.lpos+1, tmp.lcoverage);
	samclose(tmp.in);
	return 0;
}