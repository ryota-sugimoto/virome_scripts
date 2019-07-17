#!/usr/bin/env bash

[ $# == 2 ] || { echo "$(basename $0) <blastout> <out_dir>"; exit 1; }
[ -f ${1} ] || { echo "${1} not found."; exit 1; }
[ -d ${2} ] || { echo "${2} not found."; exit 1; }

blastout=${1}
out_dir=${2}

script_dir=$(cd $(dirname ${0}); pwd)
all_contigs_fasta="/home/r-sugimoto/Virome/human_gut/Mining/contigs/all/all.fasta"
bedtools="/home/r-sugimoto/tools/bedtools2/bin"
bwa="/home/r-sugimoto/tools/bwa-0.7.17/bwa"

spacers=($(cat ${blastout} | cut -f 1 | sort | uniq))

for spacer in ${spacers[@]}
do
  samtools faidx ${all_contigs_fasta} \
    $(awk -v spacer=${spacer} '$1 == spacer {print $2}' ${blastout}) \
    > ${out_dir}/tmp.fasta || exit 1
  
  ${script_dir}/split_fasta.py ${out_dir}/tmp.fasta || exit 1
  
  num_contigs=$(grep '>' ${out_dir}/tmp.fasta | wc -l)
  seq ${num_contigs} | while read n;
  do
    ${bwa} index ${out_dir}/tmp_${n}.fasta || exit 1
    
    ${bwa} mem ${out_dir}/tmp_${n}.fasta tmp_${n}_.fasta | samtools sort \
      > ${out_dir}/tmp.bam || exit 1
    
    ${bedtools}/genomeCoverageBed -ibam ${out_dir}/tmp.bam -bg \
      | ${bedtools}/mergeBed -d 500 \
      > ${out_dir}/tmp.bed || exit 1
    
    ${bedtools}/fastaFromBed -fi ${out_dir}/tmp_${n}.fasta \
                             -bed ${out_dir}/tmp.bed \
      > ${out_dir}/tmp_${n}_protospacer.fasta || exit 1
  done || exit 1
  
  cat ${out_dir}/tmp_*_protospacer.fasta \
    > ${out_dir}/${spacer}.protospacer.fasta || exit 1
  rm ${out_dir}/tmp*
done
