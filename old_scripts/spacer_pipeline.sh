#!/usr/bin/env bash

[ $# == 1 ] || { echo "$(basename ${0}) <sample_dir>"; exit 1; }
[ -d ${1} ] || { echo "${1} not exist"; exit 1; }

crisprfinder=/home/r-sugimoto/tools/CRISPRCasFinder/CRISPRCasFinder.pl
sofile=/home/r-sugimoto/tools/CRISPRCasFinder/sel392v2.so

script_dir=$(cd $(dirname ${0}); pwd)
sample_dir=$(cd ${1}; pwd)
sample_name=$(basename ${sample_dir})
fasta=${sample_dir}/sspace/standard_output.final.scaffolds.fasta
[ -f ${fasta} ] || { echo "${fasta} not found."; exit 1; }
fastq1=${sample_dir}/fastq/${sample_name}_1.clean.fastq.gz
fastq2=${sample_dir}/fastq/${sample_name}_2.clean.fastq.gz
[ -f ${fastq1} ] || { echo "${fastq1} not found."; exit 1; }
[ -f ${fastq2} ] || { echo "${fastq2} not found."; exit 1; }

min_len=10000
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

pushd $(dirname ${crisprfinder}) > /dev/null
time=$(date +%s)
mkdir -p work/${sample_name}_${time}
cd work/${sample_name}_${time}
crispr_dir=$(pwd)
cp ${processed_fasta} ./${sample_name}.fasta
cmd=(../../CRISPRCasFinder.pl 
     -cf CasFinder-2.0.2
     -def General
 #    -cas
     -i ${sample_name}.fasta
 #    -keep
     -out result
     -so ../../sel392v2.so)
${cmd[@]} &> log || exit 1
rm ${sample_name}.fasta
popd > /dev/null

mkdir ${sample_dir}/crispr
cp -r ${crispr_dir}/result/* ${sample_dir}/crispr/
crispr_dir=${sample_dir}/crispr

result_json=${crispr_dir}/result.json
[ -f ${result_json} ] || { echo "${result_json} not found"; exit 1; }

spacer_dir=${sample_dir}/spacer
mkdir ${spacer_dir}
cmd=(${script_dir}/collect_spacer.sh
     ${fastq1}
     ${fastq2}
     ${result_json}
     ${spacer_dir})
${cmd[@]} || exit 1

bed_file=${crispr_dir}/crispr.bed
cat ${crispr_dir}/TSV/Crisprs_REPORT.tsv \
  | tail -n +2 \
  | sed '/^$/d' \
  | awk 'BEGIN{OFS="\t"} 
         { if ($6 < 100) {
             begin = 0; }
           else {
             begin = $6 - 100; } 
           print $2, begin, $7+100}' \
  > ${bed_file}
masked_fasta=${processed_fasta%.fasta}.crispr_masked.fasta
cmd=(bedtools maskfasta
     -fi ${processed_fasta}
     -fo ${masked_fasta}
     -bed ${bed_file})
${cmd[@]} || exit 1

bwa index ${masked_fasta} || exit 1
