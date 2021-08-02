#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('rename_list', type=argparse.FileType('r')) 
parser.add_argument('fasta', type=argparse.FileType('r'))
args = parser.parse_args()

d = {}
for s in args.rename_list:
  l = s.strip().split()
  d[l[0]] = l[1]

from Bio import SeqIO
import sys
for record in SeqIO.parse(args.fasta, 'fasta'):
  if record.id in d:
    record.id = d[record.id]
  SeqIO.write(record, sys.stdout, 'fasta')
