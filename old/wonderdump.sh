#!/usr/bin/env bash

#
# This code is based on "Wonderdump" 
# from http://data.biostarhandbook.com/scripts/wonderdump.sh
#

echo '---------------'
[ $# == 2 ] || { echo "Usage: $(basename ${0}) <run_id> <out_dir>"; exit 1; }
[ -d ${2} ] || { echo "ERROR: ${2} not found"; exit 1; }

#TODO You must edit here
fasterq_dump="/home/r-sugimoto/tools/sratoolkit.2.9.2-ubuntu64/bin/fasterq-dump"
run_id=${1}
out_dir=$(cd ${2}; pwd)


SRA_FILE="${out_dir}/${run_id}.sra"
TMP_FILE="${out_dir}/${run_id}.tmp"

site=${run_id:0:3}
PATH1=${run_id:0:6}
PATH2=${run_id:0:10}

URL="ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/${site}/${PATH1}/${PATH2}/${run_id}.sra"
echo "URL:${URL}"
wget -nv -o /dev/stdout -O ${TMP_FILE} ${URL} || exit 1
mv $TMP_FILE $SRA_FILE

pushd ${out_dir} > /dev/null
${fasterq_dump} ${SRA_FILE} || exit 1
popd > /dev/null

pigz -p 4 -c ${out_dir}/${run_id}.sra_1.fastq > ${out_dir}/${run_id}_1.fastq.gz
pigz -p 4 -c ${out_dir}/${run_id}.sra_2.fastq > ${out_dir}/${run_id}_2.fastq.gz

rm ${out_dir}/${run_id}.sra_1.fastq
rm ${out_dir}/${run_id}.sra_2.fastq
[ -f ${out_dir}/${run_id}.sra_3.fastq ] && rm ${out_dir}/${run_id}.sra_3.fastq
rm ${SRA_FILE}
