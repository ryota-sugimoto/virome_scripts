#!/usr/bin/env bash

[ $# == 1 ] || { echo "Usage: $(basename $0) <fasta>"; exit 1; }
[ -f ${1} ] || { echo "ERROR: ${1} not found"; exit 1; }

fasta=${1}
script_dir=$(dirname ${0})
bbmerge='/home/ryota/workspace/tools/bbmap/bbmerge.sh'

broke_fasta=${fasta%.fasta}.broke.fasta
${script_dir}/break_sequence.py ${fasta} > ${broke_fasta} || exit 1

merged_fasta=${fasta%.fasta}.merged.fasta
insert=${fasta%.fasta}.insert.txt
cmd=(${bbmerge}
     t=10
     minoverlap=25
     in=${broke_fasta}
     interleaved=t
     outm=${merged_fasta}
     outinsert=${insert}
     strict=t)
${cmd[@]} || exit 1

rotated_fasta=${fasta%.fasta}.rotated.fasta
${script_dir}/rotate_sequence.py ${merged_fasta} > ${rotated_fasta} || exit 1

marked_fasta=${fasta%.fasta}.marked_circularity.fasta
${script_dir}/replace_sequence.py ${fasta} ${rotated_fasta} \
  > ${marked_fasta} || eixt 1

rm ${broke_fasta} ${merged_fasta} ${rotated_fasta} || exit 1
