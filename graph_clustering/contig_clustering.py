#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('contig_protospacers',
                    type=argparse.FileType('r'))
args = parser.parse_args()

d = {}
for s in args.contig_protospacers:
  l = s.strip().split()
  d[l[0]] = set(l[1].split(','))


def metric(s1,s2):
  return len(s1.intersection(s2))/max(len(s1),len(s2))


import networkx as nx
while d.keys():
  print('d keys', len(d.keys()))
  seed = list(d.keys())[0]
  nodes = set([])
  prev = set([seed])
  while prev:
    print('prev',len(prev))
    nodes.update(prev)
    tmp = []
    remains = d.keys() - nodes
    for node in prev:
      for key in remains:
          if len(d[node].intersection(d[key])) > 2:
            tmp.append(key)
    prev = set(tmp)
  
  for key in nodes:
    d.pop(key)
  print(len(list(nodes)))
