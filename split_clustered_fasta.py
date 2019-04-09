#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('clustered_fasta', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO

cluster = []
n = 0
for record in SeqIO.parse(args.clustered_fasta.name, 'fasta'):
  if not record.seq:
    if len(cluster) > 2:
      with open(str(n)+'.fasta', 'w') as f:
        for rec in cluster:
          SeqIO.write(rec, f, 'fasta')
    cluster = []
    n += 1
  else:
    cluster.append(record)
