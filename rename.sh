#!/usr/bin/env bash

[ $# == 1 ] || { echo "$0 <interleaved_fastq>"; exit 1; }

acc=$(gunzip -c ${1} | head -n 1 | cut -f 1 -d '.' | tr -d '@')

unpigz -p 10 -c ${1} \
  | awk -v id=${acc} \
    'BEGIN{pair=2}
     {n+=1;} \
     n%8==1{ rc+=1 } \
     n%4==1{ if (pair==2) pair=1; else pair=2; print "@"id"."rc"/"pair;} \
     n%4!=1{ print }' \
  | pigz -p 10 > ${1%.fastq.gz}.renamed.fastq.gz
