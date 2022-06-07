#!/usr/bin/env bash

[ $# == 1 ] || { echo "$(basename ${0}) <run_dir>"; exit 1; }
[ -d ${1} ] || { echo "$1 not found"; exit 1; }

run_dir=${1}
run_id=$(basename ${run_dir})

bbmap=~/tools/bbmap/bbmap.sh
bedtools=~/tools/bedtools2/bin/bedtools


dr_fasta=/home/r-sugimoto/Virome/work/human_gut_assembly/analysis/crispr/dr/all_dr.clustered.fasta
spacer_fasta=/home/r-sugimoto/Virome/work/human_gut_assembly/analysis/crispr/spacer/all_spacers.clustered.removed_short.fasta

contig_dir=${run_dir}/contig
contig_fasta=${contig_dir}/${run_id}.fasta
[ -f ${contig_fasta} ] || { echo "${contig_fasta} not found"; exit 1; }

bam=${contig_fasta%.fasta}.dr.bam
cmd=(${bbmap}
     slow=t
     k=8
     ambiguous=all
     secondary=t
     maxsites=1000
     sssr=0.8
     maxindel=1
     minid=0.93
     threads=20
     ref=${contig_fasta}
     path=${contig_dir}
     in=${dr_fasta}
     outm=${bam})
${cmd[@]} || exit 1

samtools faidx ${contig_fasta} || exit 1
genome_file=${contig_fasta%.fasta}.genome
cut -f 1,4 ${contig_fasta}.fai > ${genome_file}


mask_bed=${contig_dir}/mask.bed
${bedtools} bamtobed -i ${bam} \
  | awk 'BEGIN {FS="\t"; OFS="\t"} {split($1,a," "); $1=a[1]; print $0}' \
  | sort -k 1,1 -k2,2n -t $'\t' \
  | ${bedtools} merge -d 60 -i - \
  | ${bedtools} slop -b 100 -g ${genome_file} -i - > ${mask_bed} || exit 1


crispr_masked_fasta=${contig_fasta%.fasta}.crispr_masked.fasta
${bedtools} maskfasta -fi ${contig_fasta} -fo ${crispr_masked_fasta} -bed ${mask_bed} || exit 1
rm ${mask_bed}

bam=${contig_fasta%.fasta}.spacer.bam
cmd=(${bbmap}
     slow=t
     k=8
     ambiguous=all
     secondary=t
     maxsites=1000
     sssr=0.8
     maxindel=1
     minid=0.93
     threads=20
     ref=${crispr_masked_fasta}
     path=${contig_dir}
     in=${spacer_fasta}
     outm=${bam})
${cmd[@]} || exit 1

spacer_bed=${contig_dir}/spacer.bed
${bedtools} bamtobed -i ${bam} > ${spacer_bed} || exit 1

~/Virome/virome_scripts/analysis/spacer/filter_spacer_alignment.py \
  ${spacer_bed} \
  ${crispr_masked_fasta} \
  | sort -k 1,1 -k2,2n -t $'\t' \
  | ${bedtools} cluster -d 50000 -i - \
  | awk 'BEGIN{OFS="\t"; FS="\t"} {$1=$1"_"$7; print $0;}' \
  > ${contig_dir}/${run_id}.protospacer.bed || exit 1

rm ${spacer_bed}
rm -r ${contig_dir}/ref
chmod a-w ${contig_dir}/*
echo "${run_id} completed"
