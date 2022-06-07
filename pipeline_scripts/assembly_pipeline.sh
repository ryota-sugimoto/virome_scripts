#!/usr/bin/env bash

[ $# == 3 ] || { echo "Usage: $(basename ${0}) <run_file> <run_id> <out_dir>"; exit 1; }
[ -f ${1} ] || { echo "ERROR: ${1} not exist"; exit 1; }
[ -d ${3} ] || { echo "ERROR: ${3} not exist"; exit 1; }

script_dir=$(cd $(dirname ${0}); pwd)
run_file=$(cd $(dirname ${1}); pwd)/$(basename ${1})
run_id=${2}
out_dir=$(cd ${3}; pwd) || exit 1

#wait
#sleep $(shuf -i 1-100 -n 1)

crispr_detect="${script_dir}/crispr_detect.sh"
collect_spacer="${script_dir}/collect_spacers.sh"

#TODO You must edit here
spades=~/tools/spades-3.15.3/SPAdes-3.15.3-Linux/bin/spades.py
prodigal=~/tools/prodigal/Prodigal-2.6.3/prodigal
num_threads=20
memory_cap=100 #GiB

sample_dir=${out_dir}/${run_id}
mkdir -p ${sample_dir}/{fastq,contig,crispr,tmp} || exit 1
tmp=${sample_dir}/tmp
fastq_dir=${sample_dir}/fastq
spades_dir=${sample_dir}/spades
contig_dir=${sample_dir}/contig
crispr_dir=${sample_dir}/crispr
log=${sample_dir}/log

[ -f ${log} ] && rm ${log}

echo "Downloading ${run_id}"
pushd ${fastq_dir} > /dev/null
grep ${run_id} ${run_file} \
  | cut -f 5 | tr ';' '\n' \
  | paste - <(grep ${run_id} ${run_file} | cut -f 6 | tr ';' '\n') \
  | while read fastq_ftp fastq_md5;
    do
      wget \
        --no-verbose \
        --waitretry=30 \
        -t 20 \
        ftp://${fastq_ftp} &>> ${log} || exit 1
      echo "${fastq_md5} $(basename ${fastq_ftp})" | md5sum -c - &>> ${log} || exit 1
    done || exit 1
popd > /dev/null

[ -f ${fastq_dir}/${run_id}.fastq.gz ] \
  && rm ${fastq_dir}/${run_id}.fastq.gz

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

echo "Gene prediction ${run_id}"
${prodigal} -q -m \
            -i ${processed_fasta} \
            -a ${processed_fasta%.fasta}.proteins.fasta \
            -p meta \
            -o ${processed_fasta%.fasta}.gbk &>> ${log}

echo "Extracting spacers ${run_id}"
cp ${processed_fasta} ${crispr_dir} || exit 1
${crispr_detect} ${sample_dir} &>> ${log} || exit 1
rm ${crispr_dir}/${run_id}.fasta
${collect_spacer} ${sample_dir} &>> ${log}

makeblastdb -in ${processed_fasta} -dbtype nucl &>> ${log}
blastn \
  -outfmt 6 \
  -query ${crispr_dir}/${run_id}.crispr_dr.clustered.fasta \
  -db ${processed_fasta} \
  -qcov_hsp_perc 95 -perc_identity 95 -task 'blastn-short' 2> /dev/null \
  | cut -f 1,2 | sort | uniq \
  | sort -k 1,1 > ${crispr_dir}/${run_id}.dr_contig

rm -r ${spades_dir}/K21
rm -r ${spades_dir}/K33
rm -r ${spades_dir}/K55
pushd ${sample_dir} &> /dev/null
tar -cf - $(basename ${spades_dir}) | pigz -p ${num_threads} \
  > $(basename ${spades_dir}).tar.gz
popd &> /dev/null
rm -r ${spades_dir}
rm -r ${fastq_dir}
rm -r ${tmp}
gzip ${log}
chmod a-wx $(find ${sample_dir} -type f)
echo "Completed ${run_id}"
