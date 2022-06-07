#!/usr/bin/env bash

[ $# == 1 ] || { echo "Usage: $(basename $0) <dir>"; exit 1; }
[ -d ${1} ] || { echo "ERROR: ${1} not found."; exit 1; }

dir=${1}
id=$(basename ${dir})
fasta=${dir}/crispr/${id}.fasta

#TODO You must edit this
crispr_detect=~/tools/PROJECTS/CRISPRDetect_2.2/CRISPRDetect.pl
tmp=$(mktemp -d ~/tmp/XXXXXXXXXXXXXXXXX)

cat_fasta=${fasta%.fasta}.cat.fasta
{ echo ">${id}_allcat"; \
  egrep -v '^>' ${fasta} \
  | tr -d '\n'; } > ${cat_fasta}

[ -d /dev/shm/r-sugimoto/tmp ] || { mkdir -p /dev/shm/r-sugimoto/tmp; }

${crispr_detect} -q 0 \
                 -tmp_dir ${tmp} \
                 -f ${cat_fasta} \
                 -o ${cat_fasta%.fasta}.crispr_detect \
                 -array_quality_score_cutoff 2 || exit 1

gff=${cat_fasta%.fasta}.crispr_detect.gff
[ -f ${gff} ] || { echo 'missing gff file.'; exit 1; }


dr_n=0
> ${fasta%.fasta}.crispr_dr.fasta
grep 'repeat_region' ${gff} | cut -f 9 | while read info;
do
  ((dr_n+=1))
  dr=$(echo ${info} | tr ';' '\n' | egrep '^Note' | cut -f 2 -d '=')
  echo -e ">DR_${id}_${dr_n}\n${dr}" \
    >> ${fasta%.fasta}.crispr_dr.fasta
done 

cd-hit-est -S 1 -c 1 \
  -i ${fasta%.fasta}.crispr_dr.fasta \
  -o ${fasta%.fasta}.crispr_dr.clustered.fasta

rm ${cat_fasta}
rm -r ${tmp}
