#using c/c++ to write NGS app 
This package provides a number of small and efficient programs to perform common tasks with high throughput sequencing data in the FASTQ & BAM format. 
 
1. fastq_count: count reads for a given fastq format file with or without gzip compressed;
2. gzfastq_sample: subsample reads from a fastq format file with or without gzip compressed;
3. gzfastq_uniq: get the unique sequence and outputing to a file(s);
4. gzfastq_sort: sorting fastq files by dictionary order of read's name or sequence;
5. bam2depth: using hash to convert bam to bedGraph format,it's extremly fast with large memory consumption;
6. bam2wig: using hash to convert bam to wig format,you can also convert the wig file to bigWig format by using wig2bigWig of ucsc,it's extremly fast with large memory consumption;
7. bam_sliding_count: use a sliding window algorithm to count the read number & GC content in each window;
8. fastq_trim:cut reads to get the specified cycles;
9. pick_pair: for pair end fastq files sorted by read name, this program can pick out the reads that the read names are matched in each end.
10. Rfastqc.R: fastqc in R version.

#copying
This package is provided under a permissive MIT-style license. In particular:

Copyright (C) 2015 by Xiong Xu <xuxiong19880610@163.com> <xiongxu@me.com> 

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
