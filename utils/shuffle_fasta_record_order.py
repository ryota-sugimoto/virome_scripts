#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--seed', type=int, default=None)
parser.add_argument('fasta', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO
record_dict = SeqIO.to_dict(SeqIO.parse(args.fasta, 'fasta'))

import random
keys = list(record_dict.keys())
random.seed(args.seed)
random.shuffle(keys)

import sys
for key in keys:
  SeqIO.write(record_dict[key], sys.stdout, 'fasta')
