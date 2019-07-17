#!/usr/bin/env bash

[ $# == 1 ] || { echo "$(basename ${0}) <fasta>"; exit 1; }
[ -f ${1} ] || { echo "${1} not found"; exit 1; }

fasta=${1}

jellyfish="/home/r-sugimoto/tools/jellyfish/jellyfish-1.1.12/bin/jellyfish"

tetramers=(AAAA AAAC AAAG AAAT AACA AACC AACG AACT AAGA AAGC AAGG AAGT AATA AATC AATG AATT ACAA ACAC ACAG ACAT ACCA ACCC ACCG ACCT ACGA ACGC ACGG ACGT ACTA ACTC ACTG AGAA AGAC AGAG AGAT AGCA AGCC AGCG AGCT AGGA AGGC AGGG AGTA AGTC AGTG ATAA ATAC ATAG ATAT ATCA ATCC ATCG ATGA ATGC ATGG ATTA ATTC ATTG CAAA CAAC CAAG CACA CACC CACG CAGA CAGC CAGG CATA CATC CATG CCAA CCAC CCAG CCCA CCCC CCCG CCGA CCGC CCGG CCTA CCTC CGAA CGAC CGAG CGCA CGCC CGCG CGGA CGGC CGTA CGTC CTAA CTAC CTAG CTCA CTCC CTGA CTGC CTTA CTTC GAAA GAAC GACA GACC GAGA GAGC GATA GATC GCAA GCAC GCCA GCCC GCGA GCGC GCTA GGAA GGAC GGCA GGCC GGGA GGTA GTAA GTAC GTCA GTGA GTTA TAAA TACA TAGA TATA TCAA TCCA TCGA TGAA TGCA TTAA)

echo -e "#id\t$(echo ${tetramers[@]} | tr ' ' '\t')"
cat ${fasta} \
  | awk 'BEGIN {
           id = "";
           seq = ""; }
         /^>/ {
           if (seq) print(id"\t"seq);
           split($0,a," ");
           id = a[1];
           seq = ""; }
         !/^>/ {
           seq = seq$0;
         }
         END {
           print(id"\t"seq);
         }' \
  | while read id seq;
    do
      echo -e "${id}\n${seq}" \
        | ${jellyfish} count \
                       -m 4 \
                       -s 100M \
                       -t 4 \
                       -C \
                       -o ${fasta%.fasta}.tmp.tetramer \
                       /dev/stdin
      echo "${tetramers[@]}" | tr ' ' '\n' \
        | ${jellyfish} query \
                       ${fasta%.fasta}.tmp.tetramer_* \
        | cut -f 2 -d ' ' \
        | cat <(echo ${id} | cut -f 1 -d ' ' | tr -d '>') - \
        | tr '\n' '\t' | cat - <(echo)
      rm ${fasta%.fasta}.tmp.tetramer*
    done
