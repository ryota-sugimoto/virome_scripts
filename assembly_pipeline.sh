#!/usr/bin/env bash

[ $# == 2 ] || { echo "Usage: $(basename ${0}) <run_id> <out_dir>"; exit 1; }
[ -d ${2} ] || { echo "ERROR: ${2} not exist"; exit 1; }

#TODO You must edit here
spades="/home/ryota/workspace/tools/SPAdes/SPAdes-3.12.0-Linux/bin/spades.py"
sambamba="/home/ryota/workspace/tools/sambamba-0.6.8-linux-static"
tmp="/home/ryota/workspace/tmp"
num_threads=10
memory_cap=50

script_dir=$(cd $(dirname ${0}); pwd)

run_id=${1}
out_dir=$(cd ${2}; pwd)

mkdir -p ${out_dir}/${run_id}/{fastq,contig,result} || exit 1
fastq_dir=${out_dir}/${run_id}/fastq
contig_dir=${out_dir}/${run_id}/contig
spacer_dir=${out_dir}/${run_id}/result
log=${out_dir}/${run_id}/log

[ -f ${log} ] && rm ${log}

echo "Downloading ${run_id}"
#${script_dir}/wonderdump.sh ${run_id} ${fastq_dir} &> ${log} || exit 1
pushd ${fastq_dir} > /dev/null
/home/ryota/workspace/tools/sratoolkit.2.9.2-ubuntu64/bin/fastq-dump \
  --split-3 --gzip \
  ${run_id} &> ${log} || exit 1
popd > /dev/null

echo -e "Preprocessing ${run_id}"
fastq_1=${fastq_dir}/${run_id}_1.fastq.gz
fastq_2=${fastq_dir}/${run_id}_2.fastq.gz
${script_dir}/preprocess_pairend.sh ${fastq_1} ${fastq_2} &>> ${log} || exit 1

echo -e "Assembling ${run_id}"
assembly_cmd=(${spades}
              -t ${num_threads}
              -m ${memory_cap}
              --meta
              --only-assembler
              -1 ${fastq_dir}/${run_id}_1.ecc.fastq
              -2 ${fastq_dir}/${run_id}_2.ecc.fastq
              -o ${contig_dir})
${assembly_cmd[@]} &>> ${log} || exit 1

rm ${fastq_dir}/${run_id}_1.ecc.fastq \
   ${fastq_dir}/${run_id}_2.ecc.fastq \
   ${fastq_dir}/${run_id}_{1,2}.fastq.gz


echo -e "Extracting spacers ${run_id}"
fasta=${contig_dir}/scaffolds.fasta
min_len=1000
processed_fasta=${spacer_dir}/${run_id}.fasta
cat ${fasta} \
  | awk -v run_id=${run_id} -v min_len=${min_len} \
        'BEGIN { n = 0; }
         /^>/ { l = length(seq);
                if (n!=0 && l >= min_len) {
                  print ">CONTIG,run_id:"run_id",contig_n:"n",length:"l;
                  print seq; }
                n += 1;
                seq = ""; }
         /^[^>]/ { seq = seq""$0; }
         END { l = length(seq);
               if (length(seq) > min_len) {
                 print ">CONTIG,run_id:"run_id",contig_n:"n",length:"l;
                 print seq; }}' \
  > ${processed_fasta} || exit 1

${script_dir}/crt.sh ${processed_fasta} &>> ${log} || exit 1

${script_dir}/collect_spacer.sh \
  ${fastq_dir}/${run_id}.clean.fastq.gz \
  ${processed_fasta%.fasta}.crispr_dr.fasta \
  ${spacer_dir} &>> ${log} || exit 1

echo -e "Detecting circularity ${run_id}"
circular_fasta=${processed_fasta%.fasta}.circular.fasta
${script_dir}/circular_contigs.py \
  ${processed_fasta} \
  > ${circular_fasta} \
  2>> ${log} || exit 1

bwa index ${circular_fasta} &>> ${log} || exit 1

sam=${spacer_dir}/${run_id}.sam
alignment=(bwa mem
           -t ${num_threads}
           -p
           ${circular_fasta}
           ${fastq_dir}/${run_id}.clean.fastq.gz)
${alignment[@]} > ${sam} 2>> ${log} || exit 1

bam=${sam%.sam}.bam
convert_bam=(${sambamba} view 
             -t ${num_threads}
             -f bam
             -S
             -o ${bam}
             ${sam})
${convert_bam[@]} &>> ${log} || exit 1

sort_bam=(${sambamba} sort
          -t ${num_threads}
          --tmpdir=${tmp}
          ${bam})
${sort_bam[@]} &>> ${log} || exit 1

sorted_bam=${bam%.bam}.sorted.bam
${script_dir}/alignment_based_circularity.py \
  ${sorted_bam} \
  > ${spacer_dir}/${run_id}.pairend_circular || exit 1

rm ${sam} ${bam}

rm -r ${fastq_dir} ${contig_dir}
gzip ${log}
echo -e "Completed ${run_id}"
