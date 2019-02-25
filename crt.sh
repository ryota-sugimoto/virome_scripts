#!/usr/bin/env bash

[ $# == 1 ] || { echo "Usage: $(basename $0) <fasta>"; exit 1; }
[ -f ${1} ] || { echo "ERROR: ${0} not found."; exit 1; }

fasta=${1}

crt_jar="~/workspace/tools/crt/CRT1.2-CLI.jar"

output=${fasta%.fasta}.crtout
[ -f ${output} ] && rm ${output}
touch ${output}

tmp_crt=$(dirname ${fasta})/tmp.crtout

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
    tmp_fasta=$(dirname ${fasta})/$(echo ${id} | tr -d '>').fasta
    echo -e "${id}\n${seq}" > ${tmp_fasta}
    java -cp ~/workspace/tools/crt/CRT1.2-CLI.jar crt ${tmp_fasta} ${tmp_crt}
    echo >> ${output}
    cat ${tmp_crt} | egrep -v '^$' >> ${output}
    rm ${tmp_crt} ${tmp_fasta}
  done
