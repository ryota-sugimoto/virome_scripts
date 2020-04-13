#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--id_list',
                    type=argparse.FileType('r'))
parser.add_argument('-v', '--complement', action='store_true', default=False)
parser.add_argument('id', type=str, nargs='?')
parser.add_argument('fasta', type=argparse.FileType('r'))
args = parser.parse_args()


if args.id_list:
  id_set = set(s.strip() for s in args.id_list)
else:
  id_set = set([])

if args.id:
  id_set.update([args.id])

import sys
from Bio import SeqIO

if args.complement:
  f = lambda a,b: a not in b
else:
  f = lambda a,b: a in b

for record in SeqIO.parse(args.fasta, 'fasta'):
  if f(record.id, id_set):
    SeqIO.write(record, sys.stdout, 'fasta')
