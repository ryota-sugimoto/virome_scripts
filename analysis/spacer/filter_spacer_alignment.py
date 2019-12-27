#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('alignment_bed', type=argparse.FileType('r'))
parser.add_argument('contig_fasta', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq

contig_dict = SeqIO.index(args.contig_fasta.name, "fasta")
score_function = lambda s1,s2: pairwise2.align.globalms(s1, s2,
                                                        1, -1, -1, -1,
                                                        score_only=True)

for s in args.alignment_bed:
  l = s.split('\t')
  contig, start, end, query, score, strand = l
  start = int(start)
  end = int(end)
  query_info = dict(m.split(':') 
                    for m in query.split(' ')[1].split('=')[1].split(',')[1:])
  direct_repeat = Seq(query_info['DR_seq'])
  direct_repeat_rc = direct_repeat.reverse_complement()
  repeat_length = len(direct_repeat)
  fivedash_adjacent = contig_dict[contig][start-repeat_length: start].seq
  threedash_adjacent = contig_dict[contig][end: end+repeat_length].seq
  if len(fivedash_adjacent) == 0:
    continue
  elif len(threedash_adjacent) == 0:
    continue
  fivedash_adj_dr_score = score_function(fivedash_adjacent,
                                         direct_repeat)
  fivedash_adj_dr_rc_score = score_function(fivedash_adjacent,
                                            direct_repeat_rc)
  
  threedash_adj_dr_score = score_function(threedash_adjacent,
                                          direct_repeat)
  threedash_adj_dr_rc_score = score_function(threedash_adjacent,
                                             direct_repeat_rc)
  
  fivedash_threedash_score = score_function(fivedash_adjacent,
                                            threedash_adjacent)
  
  threshold = repeat_length*0.5
  if not any( v > threshold for v in [fivedash_adj_dr_score,
                                  fivedash_adj_dr_rc_score,
                                  threedash_adj_dr_score,
                                  threedash_adj_dr_rc_score,
                                  fivedash_threedash_score]):
    print(s.strip())
