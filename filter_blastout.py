#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
parser.add_argument('blastout', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq

record_dict = SeqIO.index(args.fasta.name, "fasta")
for s in args.blastout:
  l = s.split()
  query, subject, pidentical, length, mismath, gapopen, qstart, qend, \
  sstart, send, evalue, bitscore = l[:12]
  sstart = int(sstart)
  send = int(send)-1
  query_dict = dict([ s.split(':') for s in query.split(',')[1:] ])
  direct_repeat = Seq(query_dict['DR_seq'])
  repeat_length = len(direct_repeat)
  if sstart > send:
    sstart, send = send, sstart
  fivedash_adjacent = record_dict[subject][sstart-repeat_length: sstart].seq
  threedash_adjacent = record_dict[subject][send: send+repeat_length].seq
  score_function = lambda s1,s2: pairwise2.align.globalms(s1, s2,
                                                          1, -1, -1, -0.5,
                                                          score_only=True)
  fivedash_adj_dr_score = score_function(fivedash_adjacent, direct_repeat)
  threedash_adj_dr_score = score_function(threedash_adjacent, direct_repeat)
  fivedash_threedash_score = score_function(fivedash_adjacent,
                                            threedash_adjacent)
  
  threshold = repeat_length*0.5
  if fivedash_adj_dr_score > threshold:
    continue
  elif threedash_adj_dr_score > threshold:
    continue
  elif fivedash_threedash_score > threshold:
    continue
  else:
    print(s.strip())
