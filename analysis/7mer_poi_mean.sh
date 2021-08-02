#!/usr/bin/env bash

[ $# == 1 ] || { echo "$(basename ${0}) <fasta>"; exit 1; }
[ -f ${1} ] || { echo "${1} not found"; exit 1; }


tmpdir=$(mktemp -d ./tmp_XXXXXXXXX)
fasta=${1}
cat ${fasta} \
  | awk 'BEGIN{FS="\t"; OFS="\t"}
         /^>/ { if (id) { print id,seq; }; split($0,a," "); id=a[1]; seq=""; }
        !/^>/ { seq=seq""$0; }
          END { print id,seq; }' \
  | while read id seq;
    do
      echo -e "${id}\n${seq}" > ${tmpdir}/tmp.fasta
      jellyfish count -C -m 4 -s 100 -t 2 -o ${tmpdir}/mer ${tmpdir}/tmp.fasta
      jellyfish histo ${tmpdir}/mer \
      | awk -v id=${id} 'BEGIN{OFS="\t"; OFMT="%.15f"}
                         { sum1 += $1*$2; sum2+=$2 }
                         END{ print id,sum1,sum2,sum1/sum2 }' \
      | tr -d '>'
    done
rm -r ${tmpdir}
