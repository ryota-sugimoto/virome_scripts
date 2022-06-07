#!/usr/bin/env bash

[ $# == 1 ] || { echo 'collect_spacer.sh <dir>'; exit 1; }
[ -d ${1} ] || { echo "${1} not found."; exit 1; }

result_dir=${1}
id=$(basename ${result_dir})
fastq=${result_dir}/fastq/${id}.clean.fastq.gz
crispr_dr_fasta=${result_dir}/crispr/${id}.crispr_dr.clustered.fasta

[ -f ${fastq} ] || { echo "${fastq} not found."; exit 1; }
[ -f ${crispr_dr_fasta} ] || { echo "${crispr_dr_fasta} not found."; exit 1; }

script_dir=$(cd $(dirname ${0}); pwd)

#TODO You must edit here
bbmap=~/tools/bbtools/bbmap/
cdhit=~/tools/cd-hit-v4.6.8-2017-1208/cd-hit-est

putative_spacer_fastq=${result_dir}/crispr/tmp_putative_spacers.fastq
cmd1=(${bbmap}/bbduk.sh
       in=${fastq}
       ref=${crispr_dr_fasta}
       outm=${putative_spacer_fastq}
       interleaved=t
       k=19
       hdist=1
       rename=t)
${cmd1[@]} || exit 1

dr_num=0
cat ${crispr_dr_fasta} | paste - - | while read dr_id dr_seq;
do
  dr_num=$(( ${dr_num} + 1 ))
  tmp_dr_fasta=${result_dir}/crispr/tmp_dr.fasta
  echo ${dr_id} ${dr_seq} | tr ' ' '\n' > ${tmp_dr_fasta}
  
  filtered_fastq=${result_dir}/crispr/tmp_${dr_seq}_${dr_num}.spacer.fastq
  dr_len=$(printf ${dr_seq} | wc -c)
  [ ${dr_len} -le 25 ] && k=${dr_len} || k=25 
  cmd2=(${bbmap}/bbduk.sh
        in=${putative_spacer_fastq}
        ref=${tmp_dr_fasta}
        outm=${filtered_fastq}
        k=${k}
        hdist=1)
  ${cmd2[@]} || exit 1

  dr_masked_fastq=${filtered_fastq%.fastq}.masked.fastq
  cmd3=(${bbmap}/bbduk.sh
        in=${filtered_fastq}
        ref=${tmp_dr_fasta}
        kmask=R
        k=17
        hdist=2
        out=${dr_masked_fastq}
        mink=8)
  ${cmd3[@]} || exit 1

  spacer_fasta=${filtered_fastq%.fastq}.fasta
  ${script_dir}/extract_spacers.py \
    --run_id ${id} \
    --DR_id $(echo ${dr_id} | tr -d '>' | awk '{print $1}') \
    --DR_seq ${dr_seq} \
    ${dr_masked_fastq} > ${spacer_fasta} || exit 1
  
  ${cdhit} -M 50000 -d 1000 -s 0.9 -c 0.9 \
    -i ${spacer_fasta} \
    -o ${spacer_fasta%.fasta}.clustered.fasta || exit 1

  rm ${tmp_dr_fasta} || exit 1
done || exit 1

all_spacers_fasta=${result_dir}/crispr/${id}.all_spacers.fasta
cat ${result_dir}/crispr/*.spacer.clustered.fasta > ${all_spacers_fasta} \
  || exit 1
rm ${result_dir}/crispr/tmp_* || exit 1

${cdhit} -T 5 -M 50000 -d 1000 -sf 1 -s 0.9 -c 0.9 \
  -i ${all_spacers_fasta} \
  -o ${all_spacers_fasta%.fasta}.clustered.fasta || exit 1
