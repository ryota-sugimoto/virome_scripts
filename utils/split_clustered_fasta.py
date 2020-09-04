#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('clustered_fasta', type=argparse.FileType('r'))
parser.add_argument('-m', '--min_cluster_size', type=int, default=1)
args = parser.parse_args()

from Bio import SeqIO

cluster = []
n = 0
for record in SeqIO.parse(args.clustered_fasta.name, 'fasta'):
  if not record.seq:
    if len(cluster) > args.min_cluster_size:
      with open('{}.fasta'.format(n), 'w') as f:
        for rec in cluster:
          SeqIO.write(rec, f, 'fasta')
      n += 1
    cluster = []
  else:
    cluster.append(record)

if len(cluster) > args.min_cluster_size:
  with open('{}.fasta'.format(n), 'w') as f:
    for rec in cluster:
      SeqIO.write(rec, f, 'fasta')
