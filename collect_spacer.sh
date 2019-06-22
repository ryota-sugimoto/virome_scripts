#!/usr/bin/env bash

[ $# == 3 ] \
  || { echo 'collect_spacer.sh <fastq> <crispr_dr_fasta> <out_dir>'; exit 1; }
[ -f ${1} ] || { echo "${1} not found."; exit 1; }
[ -f ${2} ] || { echo "${2} not found."; exit 1; }
[ -d ${3} ] || { echo "${3} not found."; exit 1; }

fastq=${1}
crispr_dr_fasta=${2}
out_dir=${3}

script_dir=$(cd $(dirname ${0}); pwd)

#TODO You must edit here
bbmap="/home/r-sugimoto/tools/bbmap"
cdhit="/home/r-sugimoto/tools/cd-hit-v4.6.8-2017-1208/cd-hit-est"
#Please make sure that this is cd-hit-est, not cd-hit

clustered_dr=${out_dir}/$(basename ${crispr_dr_fasta%.fasta}.clustered.fasta)
clustered_dr_cmd=(${cdhit}
                -T 5
                -M 50000
                -s 0.9
                -c 1
                -d 1000
                -M 30000
                -i ${crispr_dr_fasta}
                -o ${clustered_dr})
${clustered_dr_cmd[@]} || exit 1

cmd1=(${bbmap}/bbduk.sh
       in=${fastq}
       ref=${clustered_dr}
       outm=${out_dir}/all.fastq
       interleaved=t
       k=19
       hdist=1
       rename=t)
${cmd1[@]} || exit 1

dr_num=0
cat ${clustered_dr} | paste - - | while read dr_id dr_seq;
do
  dr_num=$(( ${dr_num} + 1 ))
  tmp_dr_fasta=${out_dir}/tmp_dr.fasta
  echo ${dr_id} ${dr_seq} | tr ' ' '\n' > ${tmp_dr_fasta}
  
  run_id=$(printf ${dr_id} \
           | cut -f 2 -d ',' \
           | cut -f 2 -d ':')
  filtered_fastq=${out_dir}/tmp_${dr_seq}_${dr_num}.spacer.fastq
  dr_len=$(printf ${dr_seq} | wc -c)
  [ ${dr_len} -le 25 ] && k=${dr_len} || k=25 
  cmd2=(${bbmap}/bbduk.sh
        in=${out_dir}/all.fastq
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
    --run_id ${run_id} \
    --DR_seq ${dr_seq} \
    ${dr_masked_fastq} > ${spacer_fasta} || exit 1
  
  ${cdhit} -M 50000 -d 1000 -sf 1 -s 0.9 \
    -i ${spacer_fasta} \
    -o ${spacer_fasta%.fasta}.clustered.fasta || exit 1

  rm ${tmp_dr_fasta}
done || exit 1

run_id=$(head -n 1 ${clustered_dr} \
           | cut -f 2 -d ',' \
           | cut -f 2 -d ':')
all_spacers_fasta=${out_dir}/${run_id}.all_spacers.fasta
cat ${out_dir}/*.spacer.clustered.fasta > ${all_spacers_fasta}
rm ${out_dir}/tmp_*

${cdhit} -T 5 -M 50000 -d 1000 -sf 1 -s 0.9 \
  -i ${all_spacers_fasta} \
  -o ${all_spacers_fasta%.fasta}.clustered.fasta || exit 1

rm ${out_dir}/all.fastq
