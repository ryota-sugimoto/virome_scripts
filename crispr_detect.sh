#!/usr/bin/env bash

[ $# == 1 ] || { echo "Usage: $(basename $0) <fasta>"; exit 1; }
[ -f ${1} ] || { echo "ERROR: ${0} not found."; exit 1; }

fasta=${1}

crispr_detect="/home/r-sugimoto/tools/PROJECTS/CRISPRDetect_2.2/CRISPRDetect.pl"

cat_fasta=${1%.fasta}.cat.fasta
{ echo ">$(basename ${fasta})_allcat"; \
  egrep -v '^>' ${fasta} \
  | tr -d '\n'; } > ${cat_fasta}

${crispr_detect} -q 0 \
                 -f ${cat_fasta} \
                 -o ${cat_fasta%.fasta}.crispr_detect \
                 -array_quality_score_cutoff 3 || exit 1

gff=${cat_fasta%.fasta}.crispr_detect.gff
[ -f ${gff} ] || { echo 'missing gff file.'; exit 1; }


dr_n=0
grep 'repeat_region' ${gff} | cut -f 9 | while read info;
do
  ((dr_n+=1))
  dr=$(echo ${info} | tr ';' '\n' | egrep '^Note' | cut -f 2 -d '=')
  echo -e ">$(basename ${fasta})_${dr_n}\n${dr}"
done 
