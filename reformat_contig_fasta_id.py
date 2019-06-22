#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO
import sys

for record in SeqIO.parse(args.fasta, 'fasta'):
  id = record.id
  seq_type = id.split(',')[0]
  info = dict(kv.split(':') for kv in id.split(',')[1:])
  if 'contig_n' in info:
    n = info['contig_n']
    record.description = 'info={}'.format(record.id)
    record.id = '{}_{}_{}'.format(seq_type, info['run_id'], n)
    SeqIO.write(record, sys.stdout, 'fasta')
