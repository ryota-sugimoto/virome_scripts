#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta')
args = parser.parse_args()

from Bio import SeqIO
from Bio.SeqUtils import GC_skew

for record in SeqIO.parse(args.fasta, 'fasta'):
  for skew in GC_skew(record.seq, window=100):
    print(skew)
