#!/usr/bin/env bash
[ $# == 2 ] || { echo $0 '<fastq_1> <fastq_2>'; exit 1; }
[ -f ${1} ] || { echo "$1 not exist"; exit 1; }
[ -f ${2} ] || { echo "$2 not exist"; exit 1; }

#TODO You must edit here
bbmap_dir="/home/ryota/workspace/tools/bbmap"
human_ref="/home/ryota/workspace/tools/bbmap/resources/human_masked"
tmp_dir="/home/ryota/workspace/tmp"
num_threads=10
memory_cap=50

fastq_1=${1}
fastq_2=${2}

#adapter trim
trimmed_fastq=${fastq_1%_1.fastq.gz}.trimmed.fastq.gz
trim_command=(${bbmap_dir}/bbduk.sh
               t=${num_threads}
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
human_clean_fastq_1=${fastq_1%.fastq.gz}.clean.fastq
human_clean_fastq_2=${fastq_2%.fastq.gz}.clean.fastq
#human_unclean_fastq=${trimmed_fastq%.fastq.gz}.human_unclean.fastq.gz
human_decontamination_command=(${bbmap_dir}/bbmap.sh
                               interleaved=t
                               t=${num_threads}
                               path=${human_ref}
                               in=${trimmed_fastq}
                               outu1=${human_clean_fastq_1}
                               outu2=${human_clean_fastq_2}
                               #outm=${human_unclean_fastq}
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
normalized_fastq=${fastq_1%_1.fastq.gz}.normalized.fastq.gz
normalize_command=(${bbmap_dir}/bbnorm.sh
                           t=${num_threads}
                           -Xmx${memory_cap}g
                           tmpdir=${tmp_dir}
                           in1=${human_clean_fastq_1}
                           in2=${human_clean_fastq_2}
                           out=${normalized_fastq}
                           target=40
                           min=5)
${normalize_command[@]} || exit 1

#error correction
corrected_fastq=${fastq_1%_1.fastq.gz}.ecc.fastq.gz
ecc_command=(${bbmap_dir}/tadpole.sh
                     threads=${num_threads}
                     -Xmx${memory_cap}g
                     tmpdir=${tmp_dir}
                     interleaved=t
                     in=${normalized_fastq}
                     out=${corrected_fastq}
                     mode=correct)
${ecc_command[@]} || exit 1

#merge pair
merged_fastq=${fastq_1%_1.fastq.gz}.merged.fastq.gz
unmerged_fastq=${fastq_1%_1.fastq.gz}.unmerged.fastq.gz
insert_hist=${merged_fastq%.fastq.gz}.hist.txt
merge_command=(${bbmap_dir}/bbmerge.sh
               t=${num_threads}
               interleaved=t
               in=${corrected_fastq}
               out=${merged_fastq}
               outu=${unmerged_fastq}
               ihist=${insert_hist})
${merge_command[@]} || exit 1

rm ${trimmed_fastq} ${normalized_fastq} ${corrected_fastq}
