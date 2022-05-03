#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO

for record in SeqIO.parse(args.fasta, 'fasta'):
  print(record.id + "\t" + str(len(record.seq)))
