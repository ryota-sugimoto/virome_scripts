#!/usr/bin/env bash

[ $# == 4 ] \
  || { echo 'collect_spacer.sh <fastq1> <fastq2> <result_json> <out_dir>'; exit 1; }
[ -f ${1} ] || { echo "${1} not found."; exit 1; }
[ -f ${2} ] || { echo "${2} not found."; exit 1; }
[ -f ${3} ] || { echo "${3} not found."; exit 1; }
[ -d ${4} ] || { echo "${4} not found."; exit 1; }

fastq1=${1}
fastq2=${2}
result_json=${3}
out_dir=${4}

ref=${out_dir}/crispr_dr_flank.fasta
$(dirname ${0})/crispr_DR_flank.py ${result_json} > ${ref} || exit 1

bbmap=/home/r-sugimoto/tools/bbmap
cdhit=/home/r-sugimoto/tools/cd-hit-v4.6.8-2017-1208/cd-hit-est

dr_ref=${out_dir}/dr.fasta
cat ${ref} | paste - - | grep CRISPR_DR | tr '\t' '\n' > ${dr_ref} || exit 1

#collect spacers
echo 'collecting spacers...'
cmd1=(${bbmap}/bbduk.sh 
       in1=${fastq1}
       in2=${fastq2}
       ref=${dr_ref}
       outm=${out_dir}/all.fastq.gz
       interleaved=t
       k=21 
       hdist=1
       rename=t)
${cmd1[@]} || exit 1

cat ${ref} | paste - - - - - - \
  | while read dr_id dr_seq left_id left_seq right_id right_seq;
do
  dr_tmp_ref=${out_dir}/dr_tmp.fasta
  echo ${dr_id} ${dr_seq} | tr ' ' '\n' > ${dr_tmp_ref}
  
  flank_tmp_ref=${out_dir}/flank_tmp.fasta
  touch ${flank_tmp_ref}
  left_length=$(echo ${left_seq} | tr -d 'N' |  awk '{print length($1)}')
  right_length=$(echo ${right_seq} | tr -d 'N' | awk '{print length($1)}')
  [[ ${left_seq} != 'UNKNOWN' && ${left_length} -ge 31  ]] \
    && { echo ${left_id} ${left_seq} | tr ' ' '\n' > ${flank_tmp_ref}; } 
  [[ ${right_seq} != 'UNKNOWN' && ${right_length} -ge 31 ]] \
    && { echo ${right_id} ${right_seq} | tr ' ' '\n' >> ${flank_tmp_ref}; }
  
  dr_contig=$(printf ${dr_id} | cut -f 2 -d ',' | cut -f 2 -d ':')
  crispr_n=$(printf ${dr_id} | cut -f 3 -d ',' | cut -f 2 -d ':')
  id="${dr_contig}_${crispr_n}"
  seq_length=$(echo ${dr_seq} | awk '{print length($1)}')
  [[ ${seq_length} -gt 31 ]] && k=31 || k=${seq_length}
  filtered_fastq=${out_dir}/${id}.spacer.fastq.gz
  cmd2=(${bbmap}/bbduk.sh
        in=${out_dir}/all.fastq.gz
        ref=${dr_tmp_ref}
        outm=${filtered_fastq}
        k=${k}
        hdist=1)
  ${cmd2[@]} || exit 1

  dr_masked_fastq=${out_dir}/${id}.spacer.masked.fastq.gz
  cmd3=(${bbmap}/bbduk.sh
        in=${filtered_fastq}
        ref=${dr_tmp_ref}
        kmask=R
        k=17
        hdist=2
        out=${dr_masked_fastq}
        mink=8)
  ${cmd3[@]} || exit 1

  if [ -s ${flank_tmp_ref} ]
  then
    flank_masked_fastq=${out_dir}/${id}.flank_masked.fastq.gz
    cmd4=(${bbmap}/bbduk.sh
          in=${dr_masked_fastq}
          ref=${flank_tmp_ref}
          kmask=F
          k=31
          hdist=1
          out=${flank_masked_fastq}
          mink=8)
    ${cmd4[@]} || exit 1
    mv ${flank_masked_fastq} ${dr_masked_fastq}
  fi
  
  spacer_fasta=${out_dir}/${id}.spacer.fasta
  n=1
  gunzip -c ${dr_masked_fastq} | tr '\t' ' ' | paste - - - - \
    | cut -f 2 \
    | sed 's/^[ATGCNF][ATGCNF]*R/R/' \
    | sed 's/R[ATGCNF][ATGCNF]*$/R/' \
    | sed 's/RR*/R\nR/g' \
    | egrep '^R[^R][^R]*R$' \
    | tr -d 'R' \
    | egrep -v '[NF]' \
    | while read seq;
      do
        echo "$( printf ${dr_id} | sed 's/CRISPR_DR/CRISPR_spacer/' ),spacer_n:${n},DR_seq:${dr_seq}"
        echo ${seq}
        n=$(( n + 1 ))
      done > ${spacer_fasta}
  
  length_median=$(cat ${spacer_fasta} | paste - - | awk '{print length($2)}' \
                    | sort | uniq -c | awk '{print $1"\t"$2}' \
                    | sort -n -k 1,1 -r | head -n 1 | cut -f 2)
  
  echo 'spacer median:' ${length_median}
  
  filtered_spacer_fasta=${spacer_fasta%.fasta}.filtered.fasta
  cat ${spacer_fasta} | paste - - \
    | awk -v m=${length_median} 'length($2) >= 0.8*m && length($2) <= 1.2*m {print}' \
    | tr '\t' '\n' > ${filtered_spacer_fasta}
  ${cdhit} -d 1000 -sf 1 -i ${filtered_spacer_fasta} -o ${filtered_spacer_fasta%.fasta}.clustered.fasta || exit 1

  rm ${dr_tmp_ref} ${flank_tmp_ref}
done

rm ${out_dir}/all.fastq.gz ${dr_ref}
