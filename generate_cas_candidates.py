#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
parser.add_argument('tblout', type=argparse.FileType('r'))
args = parser.parse_args()

import re
from Bio import SeqIO

record_dict = SeqIO.to_dict(SeqIO.parse(args.fasta.name, 'fasta'))

import sys

reobj = re.compile('\ +')
for s in args.tblout:
  if s[0] == '#':
    continue
  l = reobj.split(s.strip())
  cas1_n = int(l[3].split('_')[-1])
  begin = 1 if cas1_n - 20 < 1 else cas1_n - 20
  end = cas1_n + 20
  for i in range(begin, end+1):
    m = l[3].split('_')
    m[-1] = str(i)
    candidate_id = "_".join(m)
    record = record_dict.get(candidate_id)
    if record:
      SeqIO.write(record, sys.stdout, "fasta")
