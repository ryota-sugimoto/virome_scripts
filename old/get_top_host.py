#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('contig_spacer', type=argparse.FileType('r'))
args = parser.parse_args()

d = {}
for s in args.contig_spacer:
  l = s.strip().split('\t')
  if len(l) != 3:
    continue
  contig = l[0]
  spacer = l[1]
  host = l[2]
  if contig in d:
    d[contig].append(host)
  else:
    d[contig] = [host]

from collections import Counter
for contig in d:
  d[contig] = Counter(d[contig])

for contig in d:
  top,c = d[contig].most_common(1)[0]
  total = sum(d[contig].values())
  print(contig, top, str(total), str(c/total), sep="\t")
