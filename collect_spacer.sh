#!/usr/bin/env bash

[ $# == 3 ] \
  || { echo 'collect_spacer.sh <fastq> <crispr_ref> <out_dir>'; exit 1; }
[ -f ${1} ] || { echo '${1} not found.'; exit 1; }
[ -f ${2} ] || { echo '${2} not found.'; exit 1; }
[ -d ${3} ] || { echo '${3} not found.'; exit 1; }

fastq=${1}
ref=${2}
out_dir=${3}

bbmap=/home/ryota/workspace/tools/bbmap
cdhit=/home/ryota/workspace/tools/cd-hit/cd-hit-v4.6.8-2017-1208/cd-hit

dr_ref=${out_dir}/dr.fasta
cat ${ref} | paste - - | grep DR | tr '\t' '\n' > ${dr_ref}

#collect spacers
echo 'collecting spacers...'
cmd1=(${bbmap}/bbduk.sh 
       in=${fastq}
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
  
  id=$(printf ${dr_id} | sed 's/^>//')
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
        echo '>'${id}'_spacer_'${n}
        echo ${seq}
        n=$(( n + 1 ))
      done > ${spacer_fasta}

  
  
  #${cdhit} -sf 1 -s 0.9 -i ${spacer_fasta} -o ${spacer_fasta%.fasta}.clustered.fasta || exit 1
    
  length_median=$(cat ${spacer_fasta} | paste - - | awk '{print length($2)}' \
                    | sort | uniq -c | awk '{print $1"\t"$2}' \
                    | sort -n -k 1,1 -r | head -n 1 | cut -f 2)
  
  echo 'spacer median: ' ${length_median}
  
  filtered_spacer_fasta=${spacer_fasta%.fasta}.filtered.fasta
  cat ${spacer_fasta} | paste - - \
    | awk -v m=${length_median} 'length($2) >= 0.8*m && length($2) <= 1.2*m {print}' \
    | tr '\t' '\n' > ${filtered_spacer_fasta}
  ${cdhit} -sf 1 -i ${filtered_spacer_fasta} -o ${filtered_spacer_fasta%.fasta}.clustered.fasta || exit 1

  rm ${dr_tmp_ref} ${flank_tmp_ref}
done

rm ${out_dir}/all.fastq.gz ${dr_ref}

#collect protospacers
echo "collecting protospacers..."
all_spacers=${out_dir}/all_spacers.fasta
cat ${out_dir}/*.filtered.clustered.fasta > ${all_spacers}

cmd=(${bbmap}/bbduk.sh
     interleaved=t
     in=${fastq}
     outm=${out_dir}/all_protospacers.fastq.gz
     ref=${all_spacers}
     k=25
     hdist=1)
${cmd[@]} || exit 1

ls ${out_dir}/*.filtered.clustered.fasta | while read spacers;
do
  putative=${spacers%.fasta}.p.fastq.gz
  cmd=(${bbmap}/bbduk.sh
       interleaved=t
       in=${out_dir}/all_protospacers.fastq.gz
       ref=${spacers}
       outm=${putative}
       k=27
       hdist=1
       rename=t)
  ${cmd[@]} || exit 1

  dr_name=$(basename ${spacers} | cut -f 1 -d '.')
  dr_fasta=${out_dir}/dr.fasta
  grep -A 1 ${dr_name} ${ref} > ${dr_fasta}

  final=${spacers%.spacer.filtered.clustered.fasta}.protospacer.fastq.gz
  cmd=(${bbmap}/bbduk.sh
       interleaved=t
       in=${putative}
       ref=${dr_fasta}
       outu=${final}
       k=17
       hdist=1)
  ${cmd[@]} || exit 1

  n_all_spacers=$(grep '>' ${spacers} | wc -l)
  n_hit_spacers=$(gunzip -c ${final} \
                  | awk 'n%4==0 {print} {n+=1}' \
                  | tr ' ' '_' \
                  | awk 'NF!=1 {print}' \
                  | cut --complement -f 1 \
                  | tr '\t' '\n' \
                  | sed 's/\=.*$//' \
                  | sort | uniq | wc -l)
  [ ${n_all_spacers} -eq 0 ] \
    && ratio=0 \
    || ratio=$(bc -l <<< "${n_hit_spacers}/${n_all_spacers}")
  id=$(basename ${final} | cut -f 1 -d '.')
  echo -e "protospacer_summary\t${id}\t${n_all_spacers}\t${n_hit_spacers}\t${ratio}"
  
  rm ${putative} ${dr_fasta}
done

rm ${all_spacers} ${out_dir}/all_protospacers.fastq.gz
