#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
args = parser.parse_args()

import re
exp = re.compile('[gcGC]')
from Bio import SeqIO
for record in SeqIO.parse(args.fasta, 'fasta'):
  length = len(record.seq)
  gc_count = len(exp.findall(str(record.seq)))
  print(record.id, gc_count/length, sep='\t')
