#!/usr/bin/env bash
[ $# == 2 ] || { echo "extract_mge.sh <cluster_file> <out_dir>"; exit 1; }
[ -f ${1} ] || { echo "${1} not found."; exit 1; }
[ -d ${2} ] || { echo "${2} not found."; exit 1; }

cluster_file=${1}
out_dir=${2}

protospacer_blastout=/home/r-sugimoto/Virome/Analysis/comprehensive/pilot/mge/all_protospacers.filtered.blastout
contigs=/home/r-sugimoto/Virome/Analysis/comprehensive/pilot/contig/all.fasta

cdhit=/home/r-sugimoto/tools/cd-hit-v4.6.8-2017-1208/cd-hit-est
nucmer=/home/r-sugimoto/tools/mummer/MUMmer3.23/nucmer
bedtools=/home/r-sugimoto/tools/bedtools2/bin/bedtools

cut -f 1,5 ${cluster_file} | while read n cluster;
do
  result_dir=${out_dir}/${n}
  mkdir -p ${result_dir}
  echo ${cluster} | tr ',' '\n' | tr '|' '\n' | egrep '^\*' | tr -d '\*' \
    | tr '|' '\n' \
  | while read keyspacer;
  do
    keyspacer_dir=${result_dir}/${keyspacer}
    mkdir -p ${keyspacer_dir}
    log=${keyspacer_dir}/log
    samtools faidx \
      ${contigs} \
      $(awk -v keyspacer=${keyspacer} '$1==keyspacer' ${protospacer_blastout} \
        | cut -f 2 | sort | uniq) \
        > ${keyspacer_dir}/protospacer_contigs.fasta \
        2> ${log}
    ${cdhit} \
      -i ${keyspacer_dir}/protospacer_contigs.fasta \
      -o ${keyspacer_dir}/protospacer_contigs.clustered.fasta \
      &>> ${log}
    ${nucmer} \
      --maxmatch \
      --nosimplify \
      --prefix=${keyspacer_dir}/seq_seq \
      ${keyspacer_dir}/protospacer_contigs.clustered.fasta \
      ${keyspacer_dir}/protospacer_contigs.clustered.fasta \
      &>> ${log}
    $(dirname ${nucmer})/show-coords \
      -r ${keyspacer_dir}/seq_seq.delta \
      > ${keyspacer_dir}/seq_seq.coords \
      2>> ${log}
    cat ${keyspacer_dir}/seq_seq.coords \
      | awk '/^=/{f=1; next;} f==1{print}' \
      | tr '|' ' ' \
      | sed 's/  */\t/g' \
      | sed 's/^\t//' \
      | awk '$8!=$9' \
      | awk 'BEGIN{OFS="\t"} {print $8,$1,$2,$0}' \
      > ${keyspacer_dir}/coords.bed \
      2>> ${log}
    awk -v keyspacer=${keyspacer} '$1==keyspacer' ${protospacer_blastout} \
      | awk 'BEGIN{OFS="\t"}
             { if ($9<$10)
               {print $2,$9,$10,$0}
               else {print $2,$10,$9,$0}}' \
      > ${keyspacer_dir}/keyspacer.bed \
      2>> ${log}
    ${bedtools} intersect \
      -a ${keyspacer_dir}/coords.bed \
      -b ${keyspacer_dir}/keyspacer.bed \
      -wa \
      | awk '$10 > 80' | cut -f 1,2,3 | sort | uniq \
      | sort -k1,1 -k2,2n \
      | ${bedtools} merge -i - \
      > ${keyspacer_dir}/mge.bed \
      2>> ${log}
    if [ -s ${keyspacer_dir}/mge.bed ]
    then 
      awk '{print $1":"$2"-"$3}' ${keyspacer_dir}/mge.bed \
        > ${keyspacer_dir}/mge.region
      samtools faidx -r ${keyspacer_dir}/mge.region ${contigs} \
        > ${keyspacer_dir}/mge.fasta \
        2>> ${log}
      ${cdhit} \
        -c 0.8 \
        -i ${keyspacer_dir}/mge.fasta \
        -o ${keyspacer_dir}/mge.clustered.fasta \
        &>> ${log}
      cat ${keyspacer_dir}/mge.clustered.fasta \
        >> ${result_dir}/${n}.fasta
    fi
  done
  if [ -s ${result_dir}/${n}.fasta ]
  then
    ${cdhit} -i ${result_dir}/${n}.fasta \
             -o ${result_dir}/${n}.clustered.fasta \
      &> ${result_dir}/log
  fi
done
