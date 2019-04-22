#!/usr/bin/env bash

[ $# == 1 ] || { echo "Usage: $(basename ${0}) <aligned_fasta>"; exit 1; }
[ -f ${1} ] || { echo "${1} not exist"; exit 1; }

hhdir="/home/r-sugimoto/tools/hh-suite-3.2.0/build/src"
unidb="/dev/shm/r-sugimoto/uniclust30_2018_08/uniclust30_2018_08"
pdb="/dev/shm/r-sugimoto/pdb70/pdb70"

fasta=${1}

hhm=${fasta%.fasta}.hhm
hhmake=(${hhdir}/hhmake -i ${fasta}
                        -o ${hhm})
${hhmake[@]} || exit 1

a3m=${fasta%.fasta}.a3m
hhblits=(${hhdir}/hhblits -cpu 15
                          -n 2
                          -d ${unidb}
                          -i ${hhm}
                          -oa3m ${a3m})
${hhblits[@]} || exit 1

hhr=${fasta%.fasta}.hhr
hhsearch=(${hhdir}/hhsearch -cpu 15
                            -d ${pdb}
                            -i ${a3m}
                            -o ${hhr})
${hhsearch[@]} || exit 1
