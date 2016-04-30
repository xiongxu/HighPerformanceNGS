//g++ -g -O3 -Wall map_kseq.cpp -o map_kseq  -lz
#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <map>
#include <iterator>
#include "klib/kseq.h"
#include "IO_stream.h"

using namespace std;

KSTREAM_INIT(gzFile, gzread, 1048576)
__KSEQ_TYPE(gzFile)
__KSEQ_READ(static)

template<typename TString>
inline bool starts_with(const TString& str, const TString& start) {
  if (start.size() > str.size()) return false;
  return str.compare(0, start.size(), start) == 0;
}

template<typename TString>
inline bool ends_with(const TString& str, const TString& end) {
  if (end.size() > str.size()) return false;
  return std::equal(end.rbegin(), end.rend(), str.rbegin());
}

map<string,kseq_t *> &loadSequences(string fileName) {
    gzFile fa=open_input_stream(fileName.c_str());
    kstream_t *ks=ks_init(fa);
    map< string, kseq_t * > *b = new map< string, kseq_t * >;
    typedef map<string, kseq_t * >::value_type kseq_ValueType ;
    int l=0;
    while (1) {
        kseq_t *seq = (kseq_t *)calloc(1, sizeof(kseq_t));
        seq->f = ks;
        if ((l = kseq_read(seq)) <0){
            free(seq);
            break;
        }
        if(b->count(string(seq->seq.s)) == 0){
            b->insert(kseq_ValueType(string(seq->seq.s),seq));
        }
    }
    cerr << b->size() << endl;
    gzclose(fa);
    return *b;
}

void print_kbtree(const map<string,kseq_t *>& b) {
    map<string, kseq_t * >::const_iterator iter;
    for (iter=b.begin(); iter!=b.end() ; ++iter) {
        kseq_t *kseq = iter->second;
        printf("%s %s\n%s\n+\n%s\n", kseq->name.s,kseq->comment.s,kseq->seq.s,kseq->qual.s);
    }
}

int main(int argc, char *argv[]) {
    map<string, kseq_t * > b = loadSequences(argv[1]);
    print_kbtree(b);
    return 0;
}

