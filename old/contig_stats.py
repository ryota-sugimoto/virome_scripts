#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO
lengthes = []
for record in SeqIO.parse(args.fasta, 'fasta'):
  lengthes.append(len(record))

lengthes.sort()

from scipy import stats
import numpy as np
import os
cs = np.cumsum(lengthes)
N50 = lengthes[(cs > cs[-1]*0.5).nonzero()[0][0]]
max_length = lengthes[-1]
print(os.path.splitext(os.path.basename(args.fasta.name))[0],
      N50, max_length, cs[-1], sep='\t')
