#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('hhr', type=argparse.FileType('r'), nargs='+')
args = parser.parse_args()

import re
regex = re.compile('  *')
query = None
for f in args.hhr:
  for s in f:
    if s[:5] == 'Query':
      if query:
        for t in hits:
          print('\t'.join((query,str(cols))+t)) 
      query = s.strip().split()[1]
      cols = int(next(f).strip().split()[1])
      hits = []
  
    if s[0] == '>':
      l = s[1:].strip().split(' ')
      subj_id = l[0]
      subj_disc = ' '.join(l[1:])
      scores = list(map(lambda s: s.replace('%','').split('=')[1],
                    regex.split(next(f).strip())))
      hits.append((subj_id,) + tuple(scores) + (subj_disc,))

  for t in hits:
    print('\t'.join((query, str(cols))+t))
