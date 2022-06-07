#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('stockholm', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys
import os

file_n = 0
records = []
for s in args.stockholm:
  if s[0] == '#' or s.strip() == '':
    continue
  if s.strip() == '//':
    for record in records:
      SeqIO.write(record, sys.stdout, 'fasta')
    file_n += 1
    records = []
    continue
  l = s.strip().split()
  records.append(SeqRecord(Seq(l[1][:-1]), id=l[0], description=''))
