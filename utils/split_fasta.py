#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO
import os
prefix = os.path.splitext(args.fasta.name)[0]

record_dict = SeqIO.index(args.fasta.name, 'fasta')
ids = list(record_dict.keys())


for n,id in enumerate(ids, 0):
  with open("{}_{}.fasta".format(prefix, n), 'w') as f:
    SeqIO.write(record_dict[id], f, 'fasta')
