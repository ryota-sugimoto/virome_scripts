#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('cluster_file', type=argparse.FileType('r'))
parser.add_argument('spacer_file', type=argparse.FileType('r'))
args = parser.parse_args()

spacers = set(s.strip() for s in args.spacer_file)
from itertools import chain
for s in args.cluster_file:
  cluster = s.split()[4].strip()
  l = set(cluster.replace('|',',').replace('*','').split(','))
  if not l.intersection(spacers):
    print(s.strip())
