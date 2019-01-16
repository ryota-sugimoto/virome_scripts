#!/usr/bin/env python3
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO
for record in SeqIO.parse(args.fasta, 'fasta'):
  half = int(len(record.seq)/2)
  print('>' + record.id + '_broke.at_' + str(half))
  print(record.seq[half:])
  print('>' + record.id + '_broke.at_' + str(half))
  print(record.seq[:half].reverse_complement())
