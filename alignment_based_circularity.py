#!/usr/bin/env python3

import sys
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('bam', type=argparse.FileType('rb'))
parser.add_argument('-m', '--min_mapq', type=int, default=55)
parser.add_argument('-o', '--out_sam',
                    type=argparse.FileType('w'),
                    default=None)
args = parser.parse_args()

import pysam
import numpy as np
bam_file = pysam.AlignmentFile(args.bam, threads=4, require_index=True)

tl = []
ql = []
num = 1000000
i = 0
for read in bam_file:
  i += 1
  if i >= num:
    break
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

if args.out_sam:
  print(bam_file.header, file=args.out_sam)

from collections import Counter
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

  begin_reads = bam_file.fetch(contig=ref, start=0, stop=edge_length,
                               multiple_iterators=True)
  
  end_reads = bam_file.fetch(contig=ref, start=length-edge_length, stop=length,
                             multiple_iterators=True) 
  begin_nexts = []
  reverse_circular_evidence = []
  chimeric_circular_evidence = []
  for read in begin_reads:
    if read.mapping_quality < args.min_mapq:
      continue
    if not read.mate_is_unmapped \
       and not read.has_tag('SA') \
       and read.reference_name == read.next_reference_name \
       and read.next_reference_start >= length-edge_length \
       and read.is_reverse and not read.mate_is_reverse:
      reverse_circular_evidence.append(read)

    elif read.has_tag('SA') \
         and not read.is_secondary \
         and read.reference_start <= ql_mean: 
      tag = read.get_tag('SA')
      tag_l = tag.split(';')
      tag_l = list(filter(None, tag_l))
      if len(tag_l) > 1:
        continue
      SA_ref, SA_pos, SA_strand, SA_cigar, SA_mapq, SA_nm = tag_l[0].split(',')
      if read.reference_name == SA_ref \
         and int(SA_mapq) > args.min_mapq \
         and int(SA_pos) >= length-ql_mean:
        chimeric_circular_evidence.append(read)
       
    elif not read.mate_is_unmapped \
         and read.reference_name != read.next_reference_name:
      begin_nexts.append((read.next_reference_name, 
                          read.next_reference_start,
                          read))

  end_nexts = []
  forward_circular_evidence = []
  for read in end_reads:
    if not read.mate_is_unmapped \
       and not read.has_tag('SA') \
       and read.reference_name == read.next_reference_name \
       and read.next_reference_start <= edge_length \
       and not read.is_reverse and read.mate_is_reverse:
      forward_circular_evidence.append(read)

    elif read.has_tag('SA') \
         and not read.is_secondary \
         and read.reference_start >= length-ql_mean: 
      tag = read.get_tag('SA')
      tag_l = tag.split(';')
      tag_l = list(filter(None, tag_l))
      if len(tag_l) > 1:
        continue
      SA_ref, SA_pos, SA_strand, SA_cigar, SA_mapq, SA_nm = tag_l[0].split(',')
      if read.reference_name == SA_ref \
         and int(SA_mapq) > args.min_mapq \
         and int(SA_pos) <= ql_mean:
        chimeric_circular_evidence.append(read)
 
    elif not read.mate_is_unmapped \
         and read.reference_name != read.next_reference_name:
      end_nexts.append((read.next_reference_name,
                        read.next_reference_start,
                        read))

  if forward_circular_evidence or reverse_circular_evidence \
     or chimeric_circular_evidence:
    cov = bam_file.count_coverage(ref)
    dop = np.sum(cov[0]+cov[1]+cov[2]+cov[3])/length
    print('\t'.join([ref, str(dop),
                     str(len(forward_circular_evidence)),
                     str(len(reverse_circular_evidence)),
                     str(len(chimeric_circular_evidence))]))
                   
  sys.stdout.flush()
  if args.out_sam:
    for read in forward_circular_evidence:
      read.query_name += '_forward_evidence'
      print(read.to_string(), file=args.out_sam)
    for read in reverse_circular_evidence:
      read.query_name += '_reverse_evidence'
      print(read.to_string(), file=args.out_sam)
    for read in chimeric_circular_evidence:
      read.query_name += '_chimeric_evidence'
      print(read.to_string(), file=args.out_sam)

  '''
  begin_next_ref_count = Counter(next_ref_name for next_ref_name,_,_ in begin_nexts)
  begin_next_refs = { ref for ref in begin_next_ref_count 
                     if begin_next_ref_count[ref] > 10 }
  end_next_ref_count = Counter(next_ref_name for next_ref_name,_,_ in end_nexts)
  end_next_refs = { ref for ref in end_next_ref_count 
                   if end_next_ref_count[ref] > 10 }
   
  shared_next_refs = begin_next_refs & end_next_refs
  for next_ref, next_pos, read in begin_nexts + end_nexts:
    if next_ref in shared_next_refs:
      print('insert_evidence', ref, '>' , next_ref, read.to_string())
  '''
