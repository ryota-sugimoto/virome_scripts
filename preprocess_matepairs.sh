#!/usr/bin/env bash
[ $# == 2 ] || { echo $0 '<fastq_1> <fastq_2>'; exit 1; }
[ -f ${1} ] || { echo "$1 not exist"; exit 1; }
[ -f ${2} ] || { echo "$2 not exist"; exit 1; }

bbmap_dir=/home/ryota/workspace/tools/bbmap #TODO you must edit here
human_ref=/home/ryota/workspace/tools/bbmap/resources/human_masked
nextclip=/home/ryota/workspace/tools/nextclip-NextClip_v1.3.1/bin/nextclip

tmp_dir=/home/ryota_workspace/tmp

fastq_1=${1}
fastq_2=${2}

#nextclip
unpigz ${fastq_1} ${fastq_2} || exit 1
fastq_1=${fastq_1%.gz}
fastq_2=${fastq_2%.gz}

n_reads=$(cat ${fastq_1} | wc -l)
n_reads=$(( ${n_reads} / 4 ))
echo '@@@@@@@@' ${n_reads}

dir=$(dirname ${fastq_1})
nextclip_out_dir=${dir}/nextclip_out
mkdir ${nextclip_out_dir}

matepair_fastq_1=${fastq_1%.fastq}.mp.fastq
matepair_fastq_2=${fastq_2%.fastq}.mp.fastq
nextclip_cmd=(${nextclip} --input_one ${fastq_1}
                          --input_two ${fastq_2}
                          --output_prefix ${nextclip_out_dir}/out
                          --number_of_reads ${n_reads})
${nextclip_cmd[@]} &> ${nextclip_out_dir}/log || exit 1

cat ${nextclip_out_dir}/out_{A,B,C}_R1.fastq > ${matepair_fastq_1} || exit 1
cat ${nextclip_out_dir}/out_{A,B,C}_R2.fastq > ${matepair_fastq_2} || exit 1
pigz ${fastq_1} ${fastq_2} ${matepair_fastq_1} ${matepair_fastq_2} || exit 1
matepair_fastq_1=${matepair_fastq_1}.gz
matepair_fastq_2=${matepair_fastq_2}.gz
rm ${nextclip_out_dir}/out_{A,B,C,D}_R{1,2}.fastq || exit 1

#adapter,quality trim
trimmed_fastq=${matepair_fastq_1%.fastq.gz}.trimmed.fastq.gz
trim_command=(${bbmap_dir}/bbduk.sh
               t=10
               in1=${matepair_fastq_1}
               in2=${matepair_fastq_2}
               out=${trimmed_fastq}
               qtrim=rl
               trimq=20
               ref=adapters,phix
               ktrim=r
               k=23
               mink=11
               hdist=1
               tpe
               tbo
               maq=20)
${trim_command[@]} || exit 1

#human decontamination
human_clean_fastq=${trimmed_fastq%.fastq.gz}.human_clean.fastq.gz
human_unclean_fastq=${trimmed_fastq%.fastq.gz}.human_unclean.fastq.gz
human_decontamination_command=(${bbmap_dir}/bbmap.sh
                               interleaved=t
                               t=10
                               path=${human_ref}
                               in=${trimmed_fastq}
                               outu=${human_clean_fastq}
                               outm=${human_unclean_fastq}
                               minratio=0.9
                               maxindel=3
                               bwr=0.16
                               bw=12
                               fast=t
                               minhits=2
                               qtrim=r
                               trimq=10
                               untrim=t
                               idtag=t
                               printunmappedcount=t
                               kfilter=25
                               maxsites=1
                               k=14)
${human_decontamination_command[@]} || exit 1

#normalize
normalized_fastq=${human_clean_fastq%.fastq.gz}.normalized.fastq.gz
normalize_command=(${bbmap_dir}/bbnorm.sh
                           t=10
                           -Xmx100g
                           tmpdir=${tmp_dir}
                           interleaved=t
                           in=${human_clean_fastq}
                           out=${normalized_fastq}
                           target=40
                           min=5)
${normalize_command[@]} || exit 1

#error correction
corrected_fastq=${normalized_fastq%.fastq.gz}.ecc.fastq.gz
ecc_command=(${bbmap_dir}/tadpole.sh
                     threads=10
                     -Xmx100g
                     tmpdir=${tmp_dir}
                     interleaved=t
                     in=${normalized_fastq}
                     out=${corrected_fastq}
                     mode=correct)
${ecc_command[@]} || exit 1
