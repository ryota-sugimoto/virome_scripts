#!/usr/bin/env bash

[ $# == 2 ] || { echo "Usage: $(basename $0) <run_file> <out_dir>"; exit 1; }
[ -f ${1} ] || { echo "ERROR: $1 not found"; exit 1; }
[ -d ${2} ] || { echo "ERROR: $2 not found"; exit 1; }

run_file=${1}
out_dir=${2}

cat ${run_file} | while read sample_id run_id Mbp;
do
  $(dirname ${0})/assembly_pipeline.sh ${run_id} ${out_dir}
done
