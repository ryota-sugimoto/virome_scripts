#!/usr/bin/env bash
[ $# == 2 ] || { echo $0 '<fastq_1> <fastq_2>'; exit 1; }
[ -f ${1} ] || { echo "$1 not exist"; exit 1; }
[ -f ${2} ] || { echo "$2 not exist"; exit 1; }

#TODO You must edit here
bbmap_dir="/home/r-sugimoto/tools/bbtools/bbmap/"
human_ref="/home/r-sugimoto/tools/bbtools/bbmap/resources/human_masked"
tmp_dir="/home/r-sugimoto/tmp"
num_threads=10
memory_cap=100

fastq_1=${1}
fastq_2=${2}

#adapter trim
trimmed_fastq=${fastq_1%_1.fastq.gz}.trimmed.fastq.gz
trim_command=(${bbmap_dir}/bbduk.sh
               -Xmx${memory_cap}g
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
human_clean_fastq=${fastq_1%_1.fastq.gz}.clean.fastq.gz
human_decontamination_command=(${bbmap_dir}/bbmap.sh
                               -Xmx${memory_cap}g
                               t=${num_threads}
                               path=${human_ref}
                               in=${trimmed_fastq}
                               outu=${human_clean_fastq}
                               interleaved=t
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

#error correction
corrected_fastq_1=${fastq_1%.fastq.gz}.ecc.fastq
corrected_fastq_2=${fastq_2%.fastq.gz}.ecc.fastq
ecc_command=(${bbmap_dir}/tadpole.sh
                     threads=${num_threads}
                     -Xmx${memory_cap}g
                     tmpdir=${tmp_dir}
                     interleaved=t
                     in=${human_clean_fastq}
                     out1=${corrected_fastq_1}
                     out2=${corrected_fastq_2}
                     mode=correct)
${ecc_command[@]} || exit 1

rm ${trimmed_fastq} 
