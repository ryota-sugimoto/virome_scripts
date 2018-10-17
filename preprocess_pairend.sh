#!/usr/bin/env bash
[ $# == 2 ] || { echo $0 '<fastq_1> <fastq_2>'; exit 1; }
[ -f ${1} ] || { echo "$1 not exist"; exit 1; }
[ -f ${2} ] || { echo "$2 not exist"; exit 1; }

bbmap_dir=/home/r-sugimoto/tools/bbmap #TODO you must edit here
human_ref=/home/r-sugimoto/tools/bbmap/resources/human_masked
tmp_dir=/home/r-sugimoto/tmp

fastq_1=${1}
fastq_2=${2}

#adapter trim
trimmed_fastq=${1%.fastq.gz}.trimmed.fastq.gz
trim_command=(${bbmap_dir}/bbduk.sh
               t=10
               in1=${1}
               in2=${2}
               out=${trimmed_fastq}
               ftm=5
               qtrim=rl
               trimq=20
               ftl=15
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

#merge pair
merged_fastq=${corrected_fastq%.fastq.gz}.merged.fastq.gz
unmerged_fastq=${corrected_fastq%.fastq.gz}.unmerged.fastq.gz
insert_hist=${merged_fastq%.fastq.gz}.hist.txt
merge_command=(${bbmap_dir}/bbmerge.sh
               t=10
               interleaved=t
               in=${corrected_fastq}
               out=${merged_fastq}
               outu=${unmerged_fastq}
               ihist=${insert_hist})
${merge_command[@]} || exit 1

#rm ${trimmed_fastq} ${human_clean_fastq} ${normalized_fastq}
