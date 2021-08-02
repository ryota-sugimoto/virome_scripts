#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-n', '--num_repeats', type=int, default=2)
parser.add_argument('fasta')
args = parser.parse_args()

import sys
from Bio import SeqIO

for record in SeqIO.parse(args.fasta, 'fasta'):
  record.seq = record.seq * args.num_repeats
  SeqIO.write(record, sys.stdout, 'fasta')
