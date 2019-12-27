#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys

parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
parser.add_argument('--seed_length', type=int, default=23)
parser.add_argument('--slide_size', type=int, default=10)
parser.add_argument('--num_slide', type=int, default=30)
parser.add_argument('--max_redundancy_coverage', type=float, default=0.5)
parser.add_argument('--max_alignment_length', type=int, default=10000)
parser.add_argument('--make_concatemer', default=False,
                    action='store_true')
parser.add_argument('--show_alignment', default=False,
                    action='store_true')
parser.add_argument('--no_trim', default=False,
                    action='store_true')
args = parser.parse_args()

for record in SeqIO.parse(args.fasta, 'fasta'):
  for i in range(args.num_slide):
    seed_begin = i*args.slide_size
    seed_end = seed_begin + args.seed_length
    seed = record.seq[seed_begin:seed_end]
    pos = record.seq.rfind(seed, seed_end)
    alignment_length = seed_begin + len(record.seq) - pos
    if pos != -1:
      if alignment_length < len(record.seq)*args.max_redundancy_coverage \
           and alignment_length < args.max_alignment_length:
        seq1 = record.seq[:alignment_length]
        seq2 = record.seq[-alignment_length:]
        if args.show_alignment:
          alignments = pairwise2.align.globalms(seq1, seq2, 1, -1, -1, -1)
          score = alignments[0][2]
        else:
          score = pairwise2.align.globalms(seq1, seq2, 1, -1, -1, -1,
                                           score_only=True)
        if score >= alignment_length*0.9:
          if not args.no_trim:
            record.seq = record.seq[:-alignment_length]
          if args.show_alignment:
            print('alignment_length='+str(alignment_length))
            print(format_alignment(*alignments[0]))
          SeqIO.write(record, sys.stdout, "fasta")
          sys.stdout.flush()
          break
