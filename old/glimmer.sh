#!/usr/bin/env bash

[ $# == 1 ] || { echo "$(basename ${0}) <fasta>"; exit 1; }
[ -f ${1} ] || { echo "${1} not found"; exit 1; }

fasta=${1}

~/tools/glimmer/glimmer3.02/bin/long-orfs \
  -o 1000 \
  -n ${fasta} \
  ${fasta%.fasta}.orf || exit 1
echo ${fasta%.fasta}.orf

~/tools/glimmer/glimmer3.02/bin/extract \
  -t \
  ${fasta} \
  ${fasta%.fasta}.orf > ${fasta%.fasta}.train || exit 1

~/tools/glimmer/glimmer3.02/bin/build-icm \
  -r \
  ${fasta%.fasta}.icm < ${fasta%.fasta}.train || exit 1

~/tools/glimmer/glimmer3.02/bin/glimmer3 \
  -o 1000 \
  -g 110 \
  -t 30 \
  ${fasta} \
  ${fasta%.fasta}.icm \
  ${fasta%.fasta}.glimmerout || exit 1

~/Virome/virome_scripts/analysis/extact_aa_from_glimmer_prediction.py \
  ${fasta} \
  ${fasta%.fasta}.glimmerout.predict \
    > ${fasta%.fasta}.proteins.fasta || exit 1
