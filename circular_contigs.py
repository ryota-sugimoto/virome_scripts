#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio import pairwise2

parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
parser.add_argument('--seed_length', type=int, default=15)
parser.add_argument('--slide_size', type=int, default=5)
parser.add_argument('--num_slide', type=int, default=10)
parser.add_argument('--max_redundancy', type=int, default=5000)
args = parser.parse_args()

for record in SeqIO.parse(args.fasta, 'fasta'):
  seed_begin = 0
  seed_end = args.seed_length
  f = False
  for i in range(args.num_slide):
    seed = record.seq[seed_begin:seed_end]
    pos = record.seq.rfind(seed, seed_end)
    l = len(record.seq) - pos
    if pos > len(record)*0.5 and l <= args.max_redundancy:
      seq1 = record.seq[:i*args.slide_size+l]
      seq2 = record.seq[pos-i*args.slide_size:]
      score = pairwise2.align.globalms(seq1, seq2, 1, -1, -1, -0.5,
                                       score_only=True)
      if score >= len(seq1)*0.9:
        f = True
        seq = record.seq[:pos-i*args.slide_size]
        print('>'+record.id+'_nonredun_{}'.format(len(seq)))
        print(seq)
        break
      else:
        break
  if not f:
    print('>'+record.id)
    print(record.seq)
    
    seed_begin += args.slide_size
    seed_end += args.slide_size
