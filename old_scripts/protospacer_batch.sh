#!/usr/bin/env bash

[ $# == 1 ] || exit 1
[ -f ${1} ] || exit 1

sample_list=${1}
spacers=~/Virome/human_gut/Mining/all_spacers/all_spacers.clustered.fasta
script=/home/r-sugimoto/Virome/virome_scripts/protospacer_contigs.sh

cat ${sample_list} | while read d;
do
  contigs=${d}/sspace/standard_output.final.scaffolds.filtered.fasta
  masked_contigs=${d}/sspace/standard_output.final.scaffolds.filtered.crispr_masked.fasta
  out_dir=${d}/protospacer
  
  mkdir -p ${out_dir}
  cmd=(${script}
       ${contigs}
       ${masked_contigs}
       ${spacers}
       ${out_dir})
  ${cmd[@]}
done
