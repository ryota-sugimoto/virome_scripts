#!/usr/bin/env python3

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
parser.add_argument('--min_match', type=int, default=25)
args = parser.parse_args()

for record in SeqIO.parse(args.fasta, 'fasta'):
  
  seed_len = args.min_match
  seed = record.seq[:seed_len]
  pos = record.seq.rfind(seed, seed_len)
  
  while pos>0:
    if pos+seed_len == len(record.seq):
      print('>'+record.id+'_ori_{}_nonredun_{}'.format(
              len(record.seq), len(record.seq[seed_len:])))
      print(record.seq[seed_len:])
      break
    seed_len += 1
    seed = record.seq[:seed_len]
    pos = record.seq.rfind(seed, seed_len)
