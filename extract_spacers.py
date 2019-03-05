#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--run_id', type=str, default='NA')
parser.add_argument('--DR_seq', type=str, default='NA')
parser.add_argument('DR_masked_fastq', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO
import re
from collections import Counter

re_obj = re.compile('R[^R]+R')

seqs = []
for record in SeqIO.parse(args.DR_masked_fastq, 'fastq'):
  seqs.extend(re_obj.findall(str(record.seq)))

seqs = list(map(lambda s: s[1:-1], seqs))
mode_length = Counter(len(s) for s in seqs).most_common(1)[0][0]
seqs = filter(lambda s:
                len(s) > 0.8*mode_length and len(s) < 1.2*mode_length \
                and 'N' not in s,
              seqs)

for i, seq in enumerate(seqs,1):
  print('>CRISPR_spacer,run_id:{},spacer_n:{},DR_seq:{}'
        .format(args.run_id, i, args.DR_seq))
  print(seq)
