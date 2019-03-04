#!/usr/bin/env bash

[ $# == 3 ] \
  || { echo 'collect_spacer.sh <fastq> <crispr_dr_fasta> <out_dir>'; exit 1; }
[ -f ${1} ] || { echo "${1} not found."; exit 1; }
[ -f ${2} ] || { echo "${2} not found."; exit 1; }
[ -d ${3} ] || { echo "${3} not found."; exit 1; }

fastq=${1}
crispr_dr_fasta=${2}
out_dir=${3}

#TODO You must edit here
bbmap=/home/ryota/workspace/tools/bbmap
cdhit=/home/ryota/workspace/tools/cd-hit/cd-hit-v4.6.8-2017-1208/cd-hit-est

clustered_dr=${out_dir}/$(basename ${crispr_dr_fasta%.fasta}.clustered.fasta)
clustered_dr_cmd=(${cdhit}
                -T 5
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
       outm=${out_dir}/all.fastq.gz
       interleaved=t
       k=19
       hdist=1
       rename=t)
${cmd1[@]} || exit 1

dr_num=0
cat ${clustered_dr} | paste - - | while read dr_id dr_seq;
do
  dr_num=$(( ${dr_num} + 1 ))
  dr_tmp_fasta=${out_dir}/dr_tmp.fasta
  echo ${dr_id} ${dr_seq} | tr ' ' '\n' > ${dr_tmp_fasta}
  
  run_id=$(printf ${dr_id} \
           | cut -f 2 -d ',' \
           | cut -f 2 -d ':')
  filtered_fastq=${out_dir}/tmp_${dr_seq}_${dr_num}.spacer.fastq.gz
  [ $(wc -c ${dr_seq}) -le 25 ] && k=$(wc -c ${dr_seq}) || k=25 
  cmd2=(${bbmap}/bbduk.sh
        in=${out_dir}/all.fastq.gz
        ref=${dr_tmp_fasta}
        outm=${filtered_fastq}
        k=${k}
        hdist=1)
  ${cmd2[@]} || exit 1

  dr_masked_fastq=${filtered_fastq%.fastq.gz}.masked.fastq.gz
  cmd3=(${bbmap}/bbduk.sh
        in=${filtered_fastq}
        ref=${dr_tmp_fasta}
        kmask=R
        k=17
        hdist=2
        out=${dr_masked_fastq}
        mink=8)
  ${cmd3[@]} || exit 1

  spacer_fasta=${filtered_fastq%.fastq.gz}.fasta
  n=1
  gunzip -c ${dr_masked_fastq} | tr '\t' ' ' | paste - - - - \
    | cut -f 2 \
    | sed 's/^[^R][^R]*R/R/' \
    | sed 's/R[^R][^R]*$/R/' \
    | sed 's/RR*/R\nR/g' \
    | egrep '^R[^R][^R]*R$' \
    | tr -d 'R' \
    | egrep -v 'N' \
    | while read seq;
      do
        echo ">CRISPR_spacer,run_id:${run_id},spacer_n:${n},DR_seq:${dr_seq}"
        echo ${seq}
        n=$(( n + 1 ))
      done > ${spacer_fasta}
  
  length_mode=$(cat ${spacer_fasta} | paste - - | awk '{print length($2)}' \
                    | sort | uniq -c | awk '{print $1"\t"$2}' \
                    | sort -n -k 1,1 -r | head -n 1 | cut -f 2)
  
  
  filtered_spacer_fasta=${spacer_fasta%.fasta}.filtered.fasta
  cat ${spacer_fasta} | paste - - \
    | awk -v m=${length_mode} 'length($2) >= 0.8*m && length($2) <= 1.2*m {print}' \
    | tr '\t' '\n' > ${filtered_spacer_fasta}

  ${cdhit} -d 1000 -sf 1 \
    -i ${filtered_spacer_fasta} \
    -o ${filtered_spacer_fasta%.fasta}.clustered.fasta || exit 1

  rm ${dr_tmp_fasta}
done

cat ${out_dir}/*.spacer.filtered.clustered.fasta > ${out_dir}/all.spacer.fasta
rm ${out_dir}/tmp_*

${cdhit} -d 1000 -sf 1 -s 0.9 \
  -i ${out_dir}/all.spacer.fasta \
  -o ${out_dir}/all.spacer.clustered.fasta || exit 1

rm ${out_dir}/all.fastq.gz
