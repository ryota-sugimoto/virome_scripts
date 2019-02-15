#!/usr/bin/env bash

[ $# == 4 ] \
  || { echo "Usage: protospacer_contigs.sh <contigs_fasta> <crispr_masked_contigs_fasta> <spacers_fasta> <out_dir>";
       exit 1; }
[ -f ${1} ] || { echo "${1} not found"; exit 1; }
[ -f ${2} ] || { echo "${2} not found"; exit 1; }
[ -f ${3} ] || { echo "${3} not found"; exit 1; }
[ -d ${4} ] || { echo "${4} not found"; exit 1; }

contigs=${1}
crispr_masked_contigs=${2}
spacers=${3}
out_dir=${4}

cmd=(makeblastdb
     -dbtype nucl
     -in ${crispr_masked_contigs})
#${cmd[@]} || exit 1

blast_out=${out_dir}/protospacer.blastout
cmd=(blastn
     -query ${spacers}
     -db ${crispr_masked_contigs}
     -outfmt 6
     -num_threads 10)
#${cmd[@]} > ${blast_out} 2> /dev/null || exit 1

protospacer_contigs=${out_dir}/protospacer_contigs_
awk '$5 <= 1 && $6 <= 1 {print $2}' ${blast_out} \
  | sort \
  | uniq -c \
  | awk '$1 > 2 {print $2}' > ${protospacer_contigs} || exit 1

protospacer_contigs_fasta=${out_dir}/protospacer_.fasta
cmd=(samtools faidx
     ${contigs}
     $(cat ${protospacer_contigs} | tr '\n' '\t'))
${cmd[@]} > ${protospacer_contigs_fasta} || exit 1
