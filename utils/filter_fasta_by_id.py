#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--id_list',
                    type=argparse.FileType('r'))
parser.add_argument('-v', '--complement', action='store_true', default=False)
parser.add_argument('-q', '--fastq', action='store_true', default=False)

parser.add_argument('id', type=str, nargs='?')
parser.add_argument('fasta', type=argparse.FileType('r'))
args = parser.parse_args()

id_set = {}
if args.id_list:
  for s in args.id_list:
    l = s.strip().split()
    if len(l) < 1:
      id_set[l[0]] = l[1]
    else:
      id_set[l[0]] = None
else:
  id_set = {}

if not args.id_list and args.id:
  id_set[args.id] = None

import sys
from Bio import SeqIO

if args.complement:
  f = lambda a,b: a not in b
else:
  f = lambda a,b: a in b

if args.fastq:
  fmt = 'fastq'
else:
  fmt = 'fasta'
for record in SeqIO.parse(args.fasta, fmt):
  if f(record.id, id_set):
    SeqIO.write(record, sys.stdout, fmt)
