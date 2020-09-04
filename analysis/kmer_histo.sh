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
      jellyfish count -m 21 -s 100M -t 10 -o ${tmpdir}/mer ${tmpdir}/tmp.fasta
      jellyfish histo ${tmpdir}/mer_0 \
      | awk -v id=${id} 'BEGIN{a[1]=0; a[2]=0; a[3]=0; OFS="\t"}
                         { a[$1] = $2 }
                         END{ print id,a[1],a[2],a[3] }' \
      | tr -d '>'
    done
rm -r ${tmpdir}
