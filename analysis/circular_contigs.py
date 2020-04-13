#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import re
import sys

parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
parser.add_argument('--seed_length', type=int, default=31)
parser.add_argument('--slide_size', type=int, default=7)
parser.add_argument('--num_slide', type=int, default=30)
parser.add_argument('--max_redundancy_coverage', type=float, default=0.7)
parser.add_argument('--max_alignment_length', type=int, default=10000)
parser.add_argument('--create_concatemer', default=False,
                    action='store_true')
parser.add_argument('--show_alignment', default=False,
                    action='store_true')
parser.add_argument('--trim_redundancy', default=False,
                    action='store_true')
parser.add_argument('--min_contig_length', type=int, default=1000)
args = parser.parse_args()

for record in SeqIO.parse(args.fasta, 'fasta'):
  found = False
  if len(record.seq) <= args.min_contig_length:
    continue
  for i in range(args.num_slide):
    seed_begin = i*args.slide_size
    seed_end = seed_begin + args.seed_length
    seed = record.seq[seed_begin:seed_end]
    for m in re.finditer(str(seed), str(record.seq[seed_end:])):
      pos = m.start() + seed_end
      alignment_length = seed_begin + len(record.seq) - pos
      if alignment_length > len(record.seq)*args.max_redundancy_coverage:
        continue
      if alignment_length > args.max_alignment_length:
        continue
      seq1 = record.seq[:alignment_length]
      seq2 = record.seq[-alignment_length:]
      if len(seq1) == 0  or len(seq1) != len(seq2):
        continue
      if args.show_alignment:
        alignments = pairwise2.align.globalms(seq1, seq2, 1, -1, -1, -1)
        score = alignments[0][2]
      else:
        score = pairwise2.align.globalms(seq1, seq2, 1, -1, -1, -1,
                                         score_only=True)
      if score < alignment_length*0.9:
        continue
      if args.trim_redundancy or args.create_concatemer:
        #record.description = record.description.split(' ')[1]
        record.id += ':0-' + str(len(record.seq)-alignment_length)
        record.seq = record.seq[:-alignment_length]
      if args.create_concatemer:
        record.id += ':concat'
        record.seq *= 3
      if args.show_alignment:
        print('alignment_length='+str(alignment_length))
        print(format_alignment(*alignments[0]))
      SeqIO.write(record, sys.stdout, "fasta")
      sys.stdout.flush()
     
      found = True
      break
    if found:
      break
