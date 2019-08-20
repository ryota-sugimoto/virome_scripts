#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('spacer_fasta', type=argparse.FileType('r'))
parser.add_argument('contig_fasta', type=argparse.FileType('r'))
parser.add_argument('blastout', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq

spacer_dict = SeqIO.index(args.spacer_fasta.name, "fasta")
contig_dict = SeqIO.index(args.contig_fasta.name, "fasta")
score_function = lambda s1,s2: pairwise2.align.globalms(s1, s2,
                                                        1, -1, -1, -0.5,
                                                        score_only=True)

for s in args.blastout:
  l = s.split()
  query, subject, pidentical, length, mismath, gapopen, qstart, qend, \
  sstart, send, evalue, bitscore = l[:12]
  sstart = int(sstart)
  send = int(send)-1
  query_dict = dict([ s.split(':') for s in \
                      spacer_dict[query].description
                                        .split(' ')[1]
                                        .split(',')[1:] ])
  direct_repeat = Seq(query_dict['DR_seq'])
  direct_repeat_rc = direct_repeat.reverse_complement()
  if sstart > send:
    sstart, send = send, sstart
  repeat_length = len(direct_repeat)
  fivedash_adjacent = contig_dict[subject][sstart-repeat_length: sstart].seq
  threedash_adjacent = contig_dict[subject][send: send+repeat_length].seq
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
  
  threshold = repeat_length*0.3
  if any( v > threshold for v in [fivedash_adj_dr_score,
                                  fivedash_adj_dr_rc_score,
                                  threedash_adj_dr_score,
                                  threedash_adj_dr_rc_score,
                                  fivedash_threedash_score]):
    continue
  else:
    print(s.strip())
