#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
parser.add_argument('id_list', type=argparse.FileType('r'))
args = parser.parse_args()

id_set = set(s.strip() for s in args.id_list)

import sys
from Bio import SeqIO
for record in SeqIO.parse(args.fasta, 'fasta'):
  if record.id in id_set:
    SeqIO.write(record, sys.stdout, 'fasta')
