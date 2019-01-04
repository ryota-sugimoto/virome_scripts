#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('bam', type=argparse.FileType('rb'))
args = parser.parse_args()

import pysam
import numpy as np
bam_file = pysam.AlignmentFile(args.bam, threads=1, require_index=True)

tl = []
ql = []
for read in bam_file:
  if not read.is_proper_pair or read.mapping_quality <= 40:
      continue
  tl.append(read.template_length)
  ql.append(read.query_length)
tl = np.array(tl)
ql = np.array(ql)
tl_mean = np.mean(np.abs(tl))
tl_std = np.std(np.abs(tl))
ql_mean = np.mean(np.abs(ql))
ql_std = np.std(np.abs(ql))

import math
for ref,length in zip(bam_file.references, bam_file.lengths):
  edge_length = 2 * ql_mean
  three_sigma = tl_mean + 3 * tl_std
  factor = three_sigma * 10
  if length < three_sigma:
    adding_edge_length = 0
  else:
    adding_edge_length = three_sigma \
                       * (1 - math.exp((three_sigma - length)/factor))
  edge_length += adding_edge_length
  begin_reads = bam_file.fetch(contig=ref, start=0, stop=edge_length)
  end_reads = bam_file.fetch(contig=ref, start=length-edge_length) 
  
  begin_mate_reads = []
  for read in begin_reads:
    if read.mate_is_unmapped:# or read.mapping_quality <= 55:
      continue
    if read.reference_name == read.next_reference_name \
       and read.next_reference_start >= length-edge_length \
       and read.is_reverse and not read.mate_is_reverse:
      print('circular_evidence_reverse' + '\t' + read.to_string())
    elif read.reference_name != read.next_reference_name:
      begin_mate_reads.append(bam_file.mate(read))

  end_mate_reads = []
  for read in end_reads:
    if read.mate_is_unmapped:
      continue
    if read.reference_name == read.next_reference_name \
       and read.next_reference_start <= edge_length \
       and not read.is_reverse and read.mate_is_reverse:
      print('circular_evidence_forward' + '\t' + read.to_string())
    elif read.reference_name != read.next_reference_name:
      end_mate_reads.append(bam_file.mate(read))

  begin_contigs = []
  for read in begin_mate_reads:
    begin_contigs.append(read.reference_name)
  begin_contigs = set(begin_contigs)
  
  end_contigs = []
  for read in end_mate_reads:
    end_contigs.append(read.reference_name)
  end_contigs = set(end_contigs)
 
  shared_contigs = begin_contigs & end_contigs
  print('num_shared_contigs', len(shared_contigs))
  for contig in shared_contigs:
    begin_shared_contig_reads = []
    for read in begin_mate_reads:
      if read.reference_name == contig:
        begin_shared_contig_reads.append(read)
    
    end_shared_contig_reads = []
    for read in end_mate_reads:
      if read.reference_name == contig:
        end_shared_contig_reads.append(read)

    print('insert_evidence_begin', ref, contig, len(begin_shared_contig_reads))
    print('insert_evidence_end', ref, contig, len(end_shared_contig_reads))
