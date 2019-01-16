#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO
for record in SeqIO.parse(args.fasta, 'fasta'):
  l = record.id.split('_')
  if 'broke.at' in l:
    index = l.index('broke.at')
    pos = int(l[index+1])
    del l[index:index+2]
    print('>' + '_'.join(l))
    print(record.seq[-pos:] + record.seq[:len(record.seq)-pos])
  else:
    print('>' + record.id)
    print(record.seq)
