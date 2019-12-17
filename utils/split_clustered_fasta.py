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
    if len(cluster) > 1:
      with open(str(rep_id)+'.fasta', 'w') as f:
        for rec in cluster:
          SeqIO.write(rec, f, 'fasta')
    rep_id = record.id
    cluster = []
    n += 1
  else:
    cluster.append(record)

if len(cluster) > 1:
  with open(str(rep_id)+'.fasta', 'w') as f:
    for rec in cluster:
      SeqIO.write(rec, f, 'fasta')
