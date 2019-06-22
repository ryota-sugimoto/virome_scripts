#!/usr/bin/env bash

[ $# == 3 ] || { echo "$(basename ${0}) <blastout> <dr_fasta> <contig_fasta>"; \
                 exit 1; }
[ -f ${1} ] || { echo "${1} not found."; exit 1; }
[ -f ${2} ] || { echo "${2} not found."; exit 1; }
[ -f ${3} ] || { echo "${3} not found."; exit 1; }

blastout=${1}
dr_fasta=${2}
contig_fasta=${3}

spacer_ids=($(cut -f 1 ${blastout} | sort | uniq))

for spacer in ${spacer_ids[@]}
do
  run_id=$(echo ${spacer} | cut -f 2 -d ',' | cut -f 2 -d ':')
  dr_seq=$(echo ${spacer} | cut -f 4 -d ',' | cut -f 2 -d ':')
  dr_ids=($(cat ${dr_fasta} | paste - - \
            | awk -v run_id=${run_id} -v dr_seq=${dr_seq} \
              '$1 ~ run_id"," && $2 == dr_seq { print $1 }'))
  contig_ids=()
  for dr_id in ${dr_ids[@]}
  do
    contig_n=$(echo ${dr_id} | cut -f 3 -d ',' | cut -f 2 -d ':')
    contig_ids+=($(printf "CONTIG_%s_%s" ${run_id} ${contig_n}))
  done
  samtools faidx ${contig_fasta} ${contig_ids[@]}
done
