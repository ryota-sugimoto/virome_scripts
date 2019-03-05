#!/usr/bin/env bash

[ $# == 1 ] || { echo "Usage: $(basename $0) <fasta>"; exit 1; }
[ -f ${1} ] || { echo "ERROR: ${0} not found."; exit 1; }

fasta=${1}

#TODO you must edit here
crt_jar=~/workspace/tools/crt/CRT1.2-CLI.jar

script_dir=$(cd $(dirname ${0}); pwd)

crt_out=${fasta%.fasta}.crtout
[ -f ${crt_out} ] && rm ${crt_out}
touch ${crt_out}

tmp_crt=$(dirname ${fasta})/tmp.crtout
tmp_fasta=$(dirname ${fasta})/tmp.fasta

cat ${fasta} \
  | awk 'BEGIN{ id="";
                seq=""; }
     /^>/ { print id"\n"seq; 
            id = $0;
            seq = ""; }
     /^[^>]/ { seq = seq""$0; }
     END { print id"\n"seq; }' \
  | tail -n +3 | paste - - \
  | while read id seq;
  do
    echo -e "${id}\n${seq}" > ${tmp_fasta}
    cmd=(java -cp ${crt_jar} crt -maxRL 70 -maxSL 70 ${tmp_fasta} ${tmp_crt})
    ${cmd[@]} || exit 1
    echo >> ${crt_out}
    cat ${tmp_crt} | egrep -v '^$' >> ${crt_out}
    rm ${tmp_crt} ${tmp_fasta}
  done || exit 1

dr_fasta=${fasta%.fasta}.crispr_dr.fasta
${script_dir}/parse_crt.py ${crt_out} > ${dr_fasta} || exit 1
