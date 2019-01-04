#!/usr/bin/env python3

import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('result_json', type=argparse.FileType('r'))
args = parser.parse_args()

d = json.load(args.result_json)

sequences = d['Sequences']

for sequence in sequences:
  if sequence['Length'] > 1000:
    id = sequence['Id']
    if sequence['Cas']:
      print('\t'.join([id] \
            + list(cas['Type'] + ':' + str(cas['Start']) \
                     for cas in sequence['Cas'])))
