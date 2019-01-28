#!/usr/bin/env bash

[ $# == 1 ] || { echo "$(basename ${0}) <sample_dir>"; exit 1; }
[ -d ${1} ] || { echo "${1} not exist"; exit 1; }

crisprfinder=/home/ryota/workspace/tools/CRISPRCasFinder/CRISPRCasFinder.pl
sofile=/home/ryota/workspace/tools/CRISPRCasFinder/sel392v2.so

script_dir=$(cd $(dirname ${0}); pwd)
sample_dir=$(cd ${1}; pwd)
sample_name=$(basename ${sample_dir})
fasta=${sample_dir}/sspace/standard_output.final.scaffolds.fasta
[ -f ${fasta} ] || { echo "fasta not found."; exit 1; }
fastq1=${sample_dir}/fastq/${sample_name}_1.clean.fastq.gz
fastq2=${sample_dir}/fastq/${sample_name}_2.clean.fastq.gz
[ -f ${fastq1} ] || { echo "fastq1 not found."; exit 1; }
[ -f ${fastq2} ] || { echo "fastq2 not found."; exit 1; }

min_len=1000
processed_fasta=${fasta%.fasta}.filtered.fasta
cat ${fasta} \
  | awk -v sample=${sample_name} -v min_len=${min_len} \
        'BEGIN { n = 0; }   
         /^>/ { l = length(seq);
                if (n!=0 && l >= min_len) { 
                  print ">"sample"_"n"_"l;
                  print seq; }
                n += 1;
                seq = ""; } 
         /^[^>]/ { seq = seq""$0; }
         END { l = length(seq);
               if (length(seq) > min_len) {
                 print ">"sample"_"n"_"l;
                 print seq; }}' \
  > ${processed_fasta}

crispr_dir=${sample_dir}/crisprfinder
mkdir ${crispr_dir}
pushd ${crispr_dir} > /dev/null
cmd=(${crisprfinder} 
     -cf CasFinder-2.0.2
     -def General
     -cas
     -i ${processed_fasta}
     -out crisprfinder.out
     -soFile ${sofile}
     -cpuM 5)
${cmd[@]} &> log || exit 1
popd > /dev/null
result_json=${crispr_dir}/crisprfinder.out/result.json
[ -f ${result_json} ] || { echo "${result_json} not found"; exit 1; }

spacer_dir=${sample_dir}/spacer
mkdir ${spacer_dir}
cmd=(${script_dir}/collect_spacer.sh
     ${fastq1}
     ${fastq2}
     ${result_json}
     ${spacer_dir})
${cmd[@]} || exit 1
