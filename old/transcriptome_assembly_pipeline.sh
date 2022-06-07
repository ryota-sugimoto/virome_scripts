#!/usr/bin/env bash

[ $# == 2 ] || { echo "Usage: $(basename ${0}) <run_id> <out_dir>"; exit 1; }
[ -d ${2} ] || { echo "ERROR: ${2} not exist"; exit 1; }

#TODO You must edit here
spades="/home/r-sugimoto/tools/SPAdes/SPAdes-3.12.0-Linux/bin/spades.py"
sambamba="/home/r-sugimoto/tools/sambamba-0.6.9-linux-static"
fasterq_dump="/home/r-sugimoto/tools/sratoolkit.2.9.2-ubuntu64/bin/fasterq-dump"
tmp="/home/r-sugimoto/tmp"
num_threads=10
memory_cap=100

script_dir=$(cd $(dirname ${0}); pwd)

run_id=${1}
out_dir=$(cd ${2}; pwd)

mkdir -p ${out_dir}/${run_id}/{fastq,contig,result} || exit 1
fastq_dir=${out_dir}/${run_id}/fastq
contig_dir=${out_dir}/${run_id}/contig
result_dir=${out_dir}/${run_id}/result
log=${out_dir}/${run_id}/log

[ -f ${log} ] && rm ${log}

echo "Downloading ${run_id}"
#${script_dir}/wonderdump.sh ${run_id} ${fastq_dir} &> ${log} || exit 1
pushd ${fastq_dir} > /dev/null
/home/r-sugimoto/tools/sratoolkit.2.9.2-ubuntu64/bin/fastq-dump \
  --split-3 -M 0 --gzip ${run_id} \
  &> ${log} || exit 1
[ -f ${fastq_dir}/${run_id}_3.fastq.gz ] \
  && rm ${fastq_dir}/${run_id}_3.fastq.gz
#${fasterq_dump} ${run_id} &> ${log} || exit 1
popd > /dev/null
#pigz -p 10 ${fastq_dir}/${run_id}_1.fastq
#pigz -p 10 ${fastq_dir}/${run_id}_2.fastq

echo -e "Preprocessing ${run_id}"
fastq_1=${fastq_dir}/${run_id}_1.fastq.gz
fastq_2=${fastq_dir}/${run_id}_2.fastq.gz
${script_dir}/preprocess_pairend.sh ${fastq_1} ${fastq_2} &>> ${log} || exit 1

echo -e "Assembling ${run_id}"
assembly_cmd=(${spades}
              -t ${num_threads}
              -m ${memory_cap}
              --rna
              -1 ${fastq_dir}/${run_id}_1.ecc.fastq
              -2 ${fastq_dir}/${run_id}_2.ecc.fastq
              -o ${contig_dir})
${assembly_cmd[@]} &>> ${log} || exit 1

rm ${fastq_dir}/${run_id}_1.ecc.fastq \
   ${fastq_dir}/${run_id}_2.ecc.fastq \
   ${fastq_dir}/${run_id}_{1,2}.fastq.gz

fasta=${contig_dir}/transcripts.fasta
min_len=0
processed_fasta=${result_dir}/${run_id}.fasta
cat ${fasta} \
  | awk -v run_id=${run_id} -v min_len=${min_len} \
        'BEGIN { n = 0; }
         /^>/ { l = length(seq);
                if (n!=0 && l >= min_len) {
                  printf(">TRANSCRIPT_%s_%i info=TRANSCRIPT,run_id:%s,contig_n:%i,length:%i,coverage:%f\n",
                         run_id, n, run_id, n, l,cov);
                  print seq; }
                n += 1;
                seq = "";
                split($0,a,"_");
                cov=a[6]; }
         /^[^>]/ { seq = seq""$0; }
         END { l = length(seq);
               if (length(seq) > min_len) {
                 printf(">TRANSCRIPT_%s_%i info=TRANSCRIPT,run_id:%s,contig_n:%i,length:%i,coverage:%f\n",
                         run_id, n, run_id, n, l, cov);
                 print seq; }}' \
  > ${processed_fasta} || exit 1


mv ${fastq_dir}/${run_id}.clean.fastq.gz ${result_dir}/
rm -r ${fastq_dir} ${contig_dir}
gzip ${log}
echo -e "Completed ${run_id}"
