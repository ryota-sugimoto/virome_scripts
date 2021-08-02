#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m', '--min_gap_propotion',
                    type=float, default=0.0)
parser.add_argument('fasta', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO
import sys


for record in SeqIO.parse(args.fasta, 'fasta'):
  if record.seq.count('N')/float(len(record.seq)) <= args.min_gap_propotion:
    SeqIO.write(record, sys.stdout, 'fasta')
