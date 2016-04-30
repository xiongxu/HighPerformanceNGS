//gcc -g -O3 -Wall kbtree_kseq.c -o kbtree_kseq -lz
#include <stdio.h>
#include <string.h>
#include "klib/kbtree.h"
#include "klib/kseq.h"
#include "IO_stream.h"

KSTREAM_INIT(gzFile, gzread, 1048576)
__KSEQ_TYPE(gzFile)
__KSEQ_READ(static)

#define seq_cmp(a, b) (a).seq.l != (b).seq.l ? (a).seq.l - (b).seq.l : strcmp((a).seq.s, (b).seq.s)
//#define seq_cmp(a, b) strcmp((a).seq.s, (b).seq.s)
KBTREE_INIT(seq, kseq_t, seq_cmp)

kbtree_t(seq) *loadSequences(const char *fileName) {
    kbtree_t(seq) *b = kb_init(seq, KB_DEFAULT_SIZE); /* KB_DEFAULT_SIZE */
    kseq_t *p;
    int l;
    gzFile fa=open_input_stream(fileName);
    kstream_t *ks=ks_init(fa);
    while (1) {
        kseq_t *kseq = (kseq_t*)calloc(1, sizeof(kseq_t));
        kseq->f = ks;
        if ((l = kseq_read(kseq)) <0){
            free(kseq);
            break;
        }
        p = kb_getp(seq, b, kseq);
        if (!p) kb_putp(seq, b, kseq);
    }
    ks_destroy(ks);
    fprintf(stderr, "%d\n", kb_size(b));
    gzclose(fa);
    return b;
}

void print_kbtree(kbtree_t(seq) *b) {
    kseq_t *p;
    kbitr_t itr;
    kb_itr_first(seq, b, &itr);
    for (; kb_itr_valid(&itr); kb_itr_next(seq, b, &itr)) {
        p = &kb_itr_key(kseq_t, &itr);
        printf("%s %s\n%s\n+\n%s\n", p->name.s,p->comment.s,p->seq.s,p->qual.s);
    }
    kb_destroy(seq, b);
}

int main(int argc, char *argv[]) {
    kbtree_t(seq) *b = loadSequences(argv[1]);
    print_kbtree(b);
    return 0;
}

