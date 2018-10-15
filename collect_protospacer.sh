#!/usr/bin/env bash

[ $# == 4 ] \
  || { echo 'collect_protospacer.sh <contigs_fasta> <crispr_dr_fasta> <spacers_fasta> <out_dir>'; \
       exit 1; }
[ -f ${1} ] || { echo '${1} not found.'; exit 1; }
[ -f ${2} ] || { echo '${2} not found.'; exit 1; }
[ -f ${3} ] || { echo '${3} not found.'; exit 1; }
[ -d ${4} ] || { echo '${4} not found.'; exit 1; }

contigs_fasta=${1}
crispr_dr_fasta=${2}
spacers_fasta=${3}
out_dir=${4}

dr_name=$(basename ${spacers_fasta} | cut -f 1 -d '.')
dr_fasta=${out_dir}/dr.fasta
grep -A 1 ${dr_name} ${crispr_dr_fasta} > ${dr_fasta}

bbmap=/home/ryota/workspace/tools/bbmap
cdhit=/home/ryota/workspace/tools/cd-hit/cd-hit-v4.6.8-2017-1208/cd-hit

#1st filter
cmd1=(${bbmap}/bbduk.sh 
       in=${contigs_fasta}
       ref=${spacers_fasta}
       outm=${out_dir}/all.fasta
       k=27
       hdist=1
       rename=t)
${cmd1[@]} || exit 1

filtered_fasta=${out_dir}/filtered.fasta
cmd2=(${bbmap}/bbduk.sh
      in=${out_dir}/all.fasta
      ref=${dr_fasta}
      outu=${filtered_fasta}
      hdist=1
      k=27)
${cmd2[@]} || exit 1

n_all_spacers=$(grep '>' ${spacers_fasta} | wc -l)
n_hit_spacers=$(grep '>' ${filtered_fasta} \
                  | cut --complement -f 1 \
                  | tr '\t' '\n' \
                  | sed 's/\=.*$//' \
                  | sort | uniq | wc -l)

echo 'Spacer summary'
echo '#all spacers:' ${n_all_spacers}
echo '#hit spacers:' ${n_hit_spacers}
echo 'ratio:' $(bc -l <<< "${n_hit_spacers}/${n_all_spacers}")

rm ${dr_fasta}
