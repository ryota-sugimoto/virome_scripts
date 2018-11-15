#!/usr/bin/env bash

#
#This code is based on "Wonderdump" 
#from http://data.biostarhandbook.com/scripts/wonderdump.sh
#

[ $# == 2 ] || { echo "Usage: $(basename ${0}) <run_id> <out_dir>"; exit 1; }
[ -d ${2} ] || { echo "ERROR: ${2} not found"; exit 1; }

fastq_dump="/home/ryota/workspace/tools/sratoolkit.2.9.2-ubuntu64/bin/fastq-dump"

run_id=${1}
out_dir=${2}

SRA_FILE="${out_dir}/${run_id}.sra"
TMP_FILE="${out_dir}/${run_id}.tmp"

PATH1=${run_id:0:6}
PATH2=${run_id:0:10}
URL="ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${PATH1}/${PATH2}/${run_id}.sra"
wget -nv -o /dev/stdout -O ${TMP_FILE} ${URL} || exit 1
mv $TMP_FILE $SRA_FILE

${fastq_dump} --split-3 --gzip ${SRA_FILE} || exit 1
rm ${SRA_FILE}
