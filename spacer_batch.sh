#!/usr/bin/env bash

[ $# == 1 ] || { echo "Usage: $(basename $0) <sample_dir_list>"; exit 1; }
[ -f ${1} ] || { echo "ERROR: ${1} not found"; exit 1; }

sample_dir_list=${1}

cat ${sample_dir_list} | while read sample_dir;
do
  $(dirname ${0})/spacer_pipeline.sh ${sample_dir} &> ${sample_dir}/spacer.log
done
