#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('begin', type=int)
parser.add_argument('end', type=int)
parser.add_argument('fasta', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO
import sys
for i, record in enumerate(SeqIO.parse(args.fasta, 'fasta')):
  if i >= args.begin and i < args.end:
    SeqIO.write(record, sys.stdout, 'fasta')
