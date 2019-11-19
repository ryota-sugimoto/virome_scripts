#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('cluster_out', type=argparse.FileType('r'))
parser.add_argument('spacers', type=argparse.FileType('r'))
args = parser.parse_args()

spacers = set(s.strip() for s in args.spacers)
from itertools import chain
for s in args.cluster_out:
  l = set(s.strip().replace('|',',').replace('*','').split(','))
  if len(l) > 3 and not l.intersection(spacers):
    print(s.strip())
