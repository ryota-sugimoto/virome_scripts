#!/usr/bin/env bash

[ $# == 2 ] || { echo "Usage: $(basename ${0}) <run_id> <out_dir>"; exit 1; }
[ -d ${2} ] || { echo "ERROR: ${2} not exist"; exit 1; }

#TODO You must edit here
spades="/home/ryota/workspace/tools/SPAdes/SPAdes-3.12.0-Linux/bin/spades.py"
sspace="/home/ryota/workspace/tools/sspace_basic-2.1.1/SSPACE_Basic.pl"
num_threads=10
memory_cap=50

script_dir=$(cd $(dirname ${0}); pwd)

run_id=${1}
out_dir=$(cd ${2}; pwd)

mkdir -p ${out_dir}/${run_id}/{fastq,spades,sspace} || exit 1
fastq_dir=${out_dir}/${run_id}/fastq
spades_dir=${out_dir}/${run_id}/spades
sspace_dir=${out_dir}/${run_id}/sspace
log=${out_dir}/${run_id}/log

echo "Downloading ${run_id}"
${script_dir}/wonderdump.sh ${run_id} ${fastq_dir} &> ${log}  || exit 1

echo -e "\nPreprocessing ${run_id}"
fastq_1=${fastq_dir}/${run_id}_1.fastq.gz
fastq_2=${fastq_dir}/${run_id}_2.fastq.gz
${script_dir}/preprocess_pairend.sh ${fastq_1} ${fastq_2} &>> ${log} || exit 1

echo -e "\nAssembling ${run_id}"
assembly_cmd=(${spades}
              -t ${num_threads}
              -m ${memory_cap}
              --meta
              --only-assembler
              --12 ${fastq_dir}/${run_id}.unmerged.fastq.gz
              --merged ${fastq_dir}/${run_id}.merged.fastq.gz
              -o ${spades_dir})
${assembly_cmd[@]} &>> ${log} || exit 1
rm -r ${spades_dir}/split_input

echo -e "\nScaffolding"
spades_fasta=${spades_dir}/scaffolds.fasta
libraries=${sspace_dir}/libraries.txt
cat << EOF > ${libraries}
PE ${fastq_dir}/${run_id}_1.clean.fastq ${fastq_dir}/${run_id}_2.clean.fastq 500 0.25 FR
EOF
pushd ${sspace_dir} > /dev/null
sspace_cmd=(${sspace}
            -x 1
            -T 10
            -l ${libraries}
            -s ${spades_fasta})
${sspace_cmd[@]} &>> ${log} || exit 1
popd > /dev/null
rm ${sspace_dir}/reads

rm ${fastq_dir}/${run_id}.merged.fastq.gz \
   ${fastq_dir}/${run_id}.merged.hist.txt \
   ${fastq_dir}/${run_id}.unmerged.fastq.gz \
   ${fastq_dir}/${run_id}_{1,2}.fastq.gz

pigz ${fastq_dir}/${run_id}_{1,2}.clean.fastq
