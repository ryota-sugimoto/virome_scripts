#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('summary_file', type=argparse.FileType('r'))
args = parser.parse_args()

length_dict = {}
for s in args.summary_file:
  l = s.strip().split()
  run_id = l[0]
  contig_length = int(l[1])
  if run_id in length_dict:
    length_dict[run_id].append(contig_length)
  else:
    length_dict[run_id] = [contig_length]

from scipy import stats
import numpy as np
for run_id in length_dict:
  lengthes = sorted(length_dict[run_id])
  var = np.std(lengthes)
  cs = np.cumsum(lengthes)
  N50 = lengthes[(cs > cs[-1]*0.5).nonzero()[0][0]]
  max_length = max(lengthes)
  print(run_id, var, N50, max_length, cs[-1], sep='\t', flush=True)
