#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('clustered_fasta', type=argparse.FileType('r'))
parser.add_argument('--no_files', action='store_true', default=False)
args = parser.parse_args()

from Bio import SeqIO
import json
json_dict = {}

ids = []
cluster = []
representative = ''
n = 0
for record in SeqIO.parse(args.clustered_fasta.name, 'fasta'):
  if not record.seq:
    if len(cluster) > 1 and not args.no_files:
      with open(str(n)+'.fasta', 'w') as f:
        for rec in cluster:
          SeqIO.write(rec, f, 'fasta')
          ids.append(rec.id)
    json_dict[representative] = ids
    representative = record.id
    ids = []
    cluster = []
    n += 1
  else:
    cluster.append(record)
    ids.append(record.id)

del json_dict['']
print(json.dumps(json_dict))
