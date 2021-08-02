#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
parser.add_argument('--split_size', '-s',
                    type=int, default=1)
parser.add_argument('--shuffle', '-r', default=False, action='store_true')
args = parser.parse_args()

from Bio import SeqIO
import os
prefix = os.path.splitext(args.fasta.name)[0]

record_dict = SeqIO.index(args.fasta.name, 'fasta')
ids = list(record_dict.keys())

from random import shuffle
if args.shuffle:
  shuffle(ids)

indxes = range(0, len(ids), args.split_size)

for n,i in enumerate(indxes, 0):
  writing_ids = ids[i: i+args.split_size]
  with open("{}.{}.fasta".format(prefix, n), 'w') as f:
    for id in writing_ids:
      SeqIO.write(record_dict[id], f, 'fasta')
