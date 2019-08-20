#!/usr/bin/env bash

[ $# == 2 ] || { echo "Usage: $(basename ${0}) <run_id> <out_dir>"; exit 1; }
[ -d ${2} ] || { echo "ERROR: ${2} not exist"; exit 1; }

#TODO You must edit here
crispr_detect="/home/r-sugimoto/Virome/virome_scripts/analysis/spacer/crispr_detect.sh"
crispr_detect="/home/r-sugimoto/Virome/virome_scripts/analysis/spacer/collect_spacer.sh"
spades="/home/r-sugimoto/tools/SPAdes/SPAdes-3.12.0-Linux/bin/spades.py"
sambamba="/home/r-sugimoto/tools/sambamba-0.6.9-linux-static"
prefetch="/home/r-sugimoto/tools/sratoolkit.2.9.2-ubuntu64/bin/prefetch"
fastq_dump="/home/r-sugimoto/tools/sratoolkit.2.9.2-ubuntu64/bin/fastq-dumph"
tmp="/home/r-sugimoto/tmp"
num_threads=10
memory_cap=100

script_dir=$(cd $(dirname ${0}); pwd)

run_id=${1}
out_dir=$(cd ${2}; pwd)

mkdir -p ${out_dir}/${run_id}/{fastq,contig,result} || exit 1
sample_dir=${out_dir}/${run_id}
fastq_dir=${sample_dir}/fastq
spades_dir=${sample_dir}/spades
contig_dir=${sample_dir}/contig
crispr_dir=${sample_dir}/crispr
log=${sample_dir}/log

[ -f ${log} ] && rm ${log}

echo "Downloading ${run_id}"
pushd ${fastq_dir} > /dev/null
${prefetch} -O ./ ${run_id} &> ${log} || exit 1
${fastq_dump} --split-3 -M 0 --gzip ${run_id}.sra &> ${log} || exit 1
popd > /dev/null

echo "Preprocessing ${run_id}"
fastq_1=${fastq_dir}/${run_id}_1.fastq.gz
fastq_2=${fastq_dir}/${run_id}_2.fastq.gz
${script_dir}/preprocess_pairend.sh ${fastq_1} ${fastq_2} &>> ${log} || exit 1

echo "Assembling ${run_id}"
assembly_cmd=(${spades}
              -t ${num_threads}
              -m ${memory_cap}
              --meta
              --only-assembler
              -1 ${fastq_dir}/${run_id}_1.ecc.fastq
              -2 ${fastq_dir}/${run_id}_2.ecc.fastq
              -o ${spades_dir})
${assembly_cmd[@]} &>> ${log} || exit 1

rm ${fastq_dir}/${run_id}_1.ecc.fastq \
   ${fastq_dir}/${run_id}_2.ecc.fastq \
   ${fastq_dir}/${run_id}_{1,2}.fastq.gz

fasta=${spades_dir}/scaffolds.fasta
min_len=1000
processed_fasta=${contig_dir}/${run_id}.fasta
cat ${fasta} \
  | awk -v run_id=${run_id} -v min_len=${min_len} \
        'BEGIN { n = 0; }
         /^>/ { l = length(seq);
                if (n!=0 && l >= min_len) {
                  printf(">CONTIG_%s_%i info=CONTIG,run_id:%s,contig_n:%i,length:%i,coverage:%f\n",
                         run_id, n, run_id, n, l,cov);
                  print seq; }
                n += 1;
                seq = "";
                split($0,a,"_");
                cov=a[6]; }
         /^[^>]/ { seq = seq""$0; }
         END { l = length(seq);
               if (length(seq) > min_len) {
                 printf(">CONTIG_%s_%i info=CONTIG,run_id:%s,contig_n:%i,length:%i,coverage:%f\n",
                         run_id, n, run_id, n, l, cov);
                 print seq; }}' \
  > ${processed_fasta} || exit 1

echo -e "Extracting spacers ${run_id}"
pushd ${crispr_dir}
cp ${processed_fasta} ./
${crispr_detect} ./${run_id}.fasta || exit 1
rm ./${run_id}.fasta
popd
${collect_spacer} ${sample_dir} &>> ${log} || exit 1

rm -r ${spades_dir}
gzip ${log}
echo "Completed ${run_id}"
