#!/usr/bin/env bash

[ $# == 2 ] || { echo "Usage: $(basename ${0}) <run_id> <out_dir>"; exit 1; }
[ -d ${2} ] || { echo "ERROR: ${2} not exist"; exit 1; }

#TODO you must edit here
spades="/home/ryota/workspace/tools/SPAdes/SPAdes-3.12.0-Linux/bin/spades.py"
sspace="/home/ryota/workspace/tools/sspace_basic-2.1.1/SSPACE_Basic.pl"

script_dir=$(cd $(dirname ${0}); pwd)

run_id=${1}
out_dir=$(cd ${2}; pwd)

echo ${out_dir}
mkdir -p ${out_dir}/${run_id}/{fastq,spades,sspace} || exit 1
fastq_dir=${out_dir}/${run_id}/fastq
spades_dir=${out_dir}/${run_id}/spades
sspace_dir=${out_dir}/${run_id}/sspace

echo "Downloading ${run_id}"
${script_dir}/wonderdump.sh ${run_id} ${fastq_dir} || exit 1

echo -e "\nPreprocessing ${run_id}"
fastq_1=${fastq_dir}/${run_id}_1.fastq.gz
fastq_2=${fastq_dir}/${run_id}_2.fastq.gz
${script_dir}/preprocess_pairend.sh ${fastq_1} ${fastq_2} || exit 1

echo -e "\nAssembling ${run_id}"
assembly_cmd=(${spades}
              -t 10
              -m 50
              --meta
              --only-assembler
              --12 ${fastq_dir}/${run_id}.unmerged.fastq.gz
              --merged ${fastq_dir}/${run_id}.merged.fastq.gz
              -o ${spades_dir})
${assembly_cmd[@]} || exit 1

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
            -s ${spade_fasta})
${sspace_cmd[@]} || exit 1
popd > /dev/null
