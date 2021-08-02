#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO
import sys

for record in SeqIO.parse(args.fasta, 'fasta'):
  record.seq = record.seq.ungap()
  SeqIO.write(record, sys.stdout, 'fasta')
