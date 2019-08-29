#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('contig_protospacers',
                    type=argparse.FileType('r'))
args = parser.parse_args()

d = {}
for s in args.contig_protospacers:
  l = s.strip().split()
  d[l[0]] = l[1].split(',')

print(d)
