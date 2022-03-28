#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('DR_masked_fastq', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO
import re
from collections import Counter

re_obj = re.compile('R[^R]+R')

seqs = []
for record in SeqIO.parse(args.DR_masked_fastq, 'fastq'):
  seqs.extend(re_obj.findall(str(record.seq)))


if seqs:
  seqs = list(map(lambda s: s[1:-1], seqs))

for i, seq in enumerate(seqs,1):
  print('>SPACER_{}'.format(i))
  print(seq)
