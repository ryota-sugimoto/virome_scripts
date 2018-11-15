#!/usr/bin/env bash

[ $# == 3 ] || { echo "Usage: $(basename $0) gzipped_fastq fq1 fq2"; exit 1; }
[ -f ${1} ] || { echo "ERROR: $1 not found"; exit 1; }

unpigz -c ${1} \
  | paste - - - - - - - - \
  | tee >(cut -f 1-4 | tr '\t' '\n' > ${2}) \
  | cut -f 5-8 \
  | tr '\t' '\n' > ${3}
