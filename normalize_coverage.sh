#!/usr/bin/env bash
[ $# == 2 ] || { echo $0 '<fastq_1> <fastq_2>'; exit 1; }
[ -f ${1} ] || { echo "$1 not exist"; exit 1; }
[ -f ${2} ] || { echo "$2 not exist"; exit 1; }

bbmap_dir=/home/ryota/workspace/tools/bbmap #TODO you must edit here
tmp_dir=/home/ryota/workspace/tmp

fastq_1=${1}
fastq_2=${2}

#adapter trim
trimmed_fastq=${1%.fastq.gz}.trimmed.fastq
trim_command=(${bbmap_dir}/bbduk.sh
               t=10
               in1=${1}
               in2=${2}
               out=${trimmed_fastq}
               ftm=5
               ftl=10
               ref=adapters
               ktrim=r
               k=23
               mink=11
               hdist=1
               tpe
               tbo
               maq=10)
${trim_command[@]} || exit 1

#human decontamination
path="/home/ryota/workspace/db/human_masked"
human_clean_fastq=${trimmed_fastq%.fastq}.human_clean.fastq
human_unclean_fastq=${trimmed_fastq%.fastq}.human_unclean.fastq
human_decontamination_command=(${bbmap_dir}/bbmap.sh
                               interleaved=t
                               path=${path}
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

#phix decontamination
phix_clean_fastq=${human_clean_fastq%.fastq}.phix_clean.fastq
phix_unclean_fastq=${human_clean_fastq%.fastq}.phix_unclean.fastq
phix_decontamination_command=(${bbmap_dir}/bbduk.sh
               t=10
               interleaved=t
               in=${human_clean_fastq}
               out=${phix_clean_fastq}
               outm=${phix_unclean_fastq}
               ref=phix
               k=31)
${phix_decontamination_command[@]} || exit 1

#normalize
normalized_fastq=${phix_clean_fastq%.fastq}.normalized.fastq
normalize_command=(${bbmap_dir}/bbnorm.sh
                           t=10
                           -Xmx100g
                           tmpdir=${tmp_dir}
                           interleaved=t
                           in=${phix_clean_fastq}
                           out=${normalized_fastq}
                           target=100
                           min=5)
${normalize_command[@]} || exit 1

#merge pair
merged_fastq=${normalized_fastq%.fastq}.merged.fastq.gz
unmerged_fastq=${normalized_fastq%.fastq}.unmerged.fastq.gz
insert_hist=${merged_fastq%.fastq}.hist.txt
merge_command=(${bbmap_dir}/bbmerge.sh
               t=10
               interleaved=t
               in=${normalized_fastq}
               out=${merged_fastq}
               outu=${unmerged_fastq}
               ihist=${insert_hist})
${merge_command[@]} || exit 1

#rm ${trimmed_fastq} ${human_clean_fastq} ${normalized_fastq}
