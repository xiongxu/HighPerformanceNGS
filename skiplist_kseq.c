// gcc -g -O3 -Wall -Wno-unused-function skiplist_kseq.c  -o skiplist_kseq -I./skiplist -lz
#include <stdio.h>
#include <string.h>
#include "skiplist_with_rank_void.h"
#include "klib/kseq.h"
#include "IO_stream.h"

KSTREAM_INIT(gzFile, gzread, 1048576)
__KSEQ_TYPE(gzFile)
__KSEQ_READ(static)

int sklCompareKey(const void *key1,const void *key2){
    return strcmp((const char *)key1,(const char *)key2);
}

void valDestructor(const void *obj){
    kseq_t *kseq = (kseq_t*)obj;
    free(kseq->name.s);
    free(kseq->comment.s);
    free(kseq->seq.s);
    free(kseq->qual.s);
    free(kseq);
}

struct skiplist *loadSequences(const char *fileName) {
    struct skiplist *skl = skiplist_new(); /* KB_DEFAULT_SIZE */
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
        struct skipnode *node=skiplist_insert(skl,kseq->seq.s,kseq,NULL,NULL,sklCompareKey);
        if (node==NULL) {
            fprintf(stderr,"err in skiplist_insert %s\n",kseq->seq.s);
            exit(EXIT_FAILURE);
        }
    }
    ks_destroy(ks);
    fprintf(stderr, "%d\n", skl->count);
    gzclose(fa);
    return skl;
}

static unsigned int print_in_rank(struct skiplist *list, unsigned int start, unsigned int stop,void (*keyDestructor)(const void *key), void (*valDestructor)(const void *obj)) {
    int i;
    unsigned int removed = 0, traversed = 0;
    struct sk_link *pos, *n, *end; //, *update[MAX_LEVEL];
    if (start <= 0 || stop <= 0 || start > list->count)
            return 0;
    i = list->level - 1;
    pos = &list->head[i];
    end = &list->head[i];
    for (; i >= 0; i--) {
        pos = pos->next;
        skiplist_foreach_forward_safe(pos, n, end) {
            struct skipnode *node = list_entry(pos, struct skipnode, link[i]);
            int current_span = traversed + node->link[i].span;
            if (current_span > stop) {
                end = &node->link[i];
                break;
            } else if (current_span >= start && i==0) {
                kseq_t *ks=(kseq_t *)node->value;
                fprintf(stderr,"%s %s %d %d %d %u\n%s\n+\n%s\n", ks->name.s,ks->comment.s,i,traversed,node->link[i].span,skiplist_key_rank(list,node->key,sklCompareKey),ks->seq.s,ks->qual.s);
                removed++;
                //continue;
            }else {
                end = &node->link[i];
                break;
            }
            traversed += node->link[i].span;
        }
        //update[i] = end;
        pos = end->prev;
        pos--;
        end--;
    }
    fprintf(stderr,"removed %d nodes\n",removed);
    return removed;
}

void skiplist_dump_kseq(struct skiplist *list) {
    struct sk_link *n;
    struct sk_link *pos = list->head[0].next;
    skiplist_foreach_forward_safe(pos, n, &list->head[0]) {
            struct skipnode *node = list_entry(pos, struct skipnode, link[0]);
            kseq_t *ks=(kseq_t *)node->value;
            fprintf(stdout,"%s %s %u\n%s\n+\n%s\n", ks->name.s,ks->comment.s,skiplist_key_rank(list,node->key,sklCompareKey),ks->seq.s,ks->qual.s);
    }
}

int main(int argc, char *argv[]) {
    struct skiplist *list = loadSequences(argv[1]);
    // skiplist_dump_kseq(list);
    int i=0;
    for (i=1;i<=20;++i){
        struct skipnode *node = skiplist_search_by_rank(list,i);
        kseq_t *ks=(kseq_t *)node->value;
        fprintf(stdout,"%s %s %u\n%s\n+\n%s\n", ks->name.s,ks->comment.s,skiplist_key_rank(list,node->key,sklCompareKey),ks->seq.s,ks->qual.s);
    }
    // remove_in_rank(list,1,20,NULL,valDestructor);
    //print_in_rank(list,1,20,NULL,valDestructor);
    skiplist_delete(list,NULL,valDestructor);
    return 0;
}

