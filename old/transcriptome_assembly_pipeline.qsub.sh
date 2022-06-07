#!/usr/bin/env bash
#$ -S /bin/bash
#$ -N assemble_gut_transcriptome
#$ -o /home/r-sugimoto/Virome/Mining/transcriptome/human_gut/logs
#$ -e /home/r-sugimoto/Virome/Mining/transcriptome/human_gut/logs

list=/home/r-sugimoto/Virome/Mining/transcriptome/human_gut/runs
script=/home/r-sugimoto/Virome/virome_scripts/assembly/transcriptome_assembly_pipeline.sh
script_dir=$(dirname ${script})
out_dir=/home/r-sugimoto/Virome/Mining/transcriptome/human_gut

n=$(( ${SGE_TASK_ID} - 1 ))
batch_size=20

cat ${list} \
  | tail -n +$(( ${n} * ${batch_size} + 1 )) \
  | head -n ${batch_size} \
  > ${out_dir}/run_${n}

cat ${out_dir}/run_${n} | while read id;
do
  ${script} ${id} ${out_dir} || { echo "${id} failed" 1>&2; }
done
