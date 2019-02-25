#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('crt_out', type=argparse.FileType('r'))
args = parser.parse_args()

import re
from Bio.Seq import Seq
f = False
for s in args.crt_out:
  s = re.sub(' +', ' ', s.strip())
  l = s.split()
  if not s:
    continue
  if l[0] == 'ORGANISM:':
    contig = l[1]
  elif l[0] == 'CRSIPR':
    n_crispr = l[1]
    crispr_range = (l[3],l[5])
  elif s[0] == '-':
    if not f:
      DRs = []
      f = True
    if f:
      consensus_dr = consensus(DRs)
      print(">CRISPR_DR,contig:{},position:{}-{}".format((contig,)+crispr_range)
      print(consensus_dr)
      f = False
  elif f:
    DRs.append(l[1])
