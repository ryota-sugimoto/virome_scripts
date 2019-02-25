#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta1', type=argparse.FileType('r'))
parser.add_argument('fasta2', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO
seq_id = []
reads = {}
for record in SeqIO.parse(args.fasta2, 'fasta'):
  seq_id.append(record.id)
  reads[record.id] = record.seq
seq_id = set(seq_id)

args.fasta2.seek(0)
for record in SeqIO.parse(args.fasta1, 'fasta'):
  if record.id in seq_id and len(reads[record.id])/len(record.seq) > 0.5:
    print('>' + record.id + '_circularity.contig_1_size_' \
          + str(len(reads[record.id])))
    print(reads[record.id])
  else:
    print('>' + record.id + '_circularity.contig_0')
    print(record.seq)
