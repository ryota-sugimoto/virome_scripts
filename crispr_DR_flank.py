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
    crisprs = sequence['Crisprs']
    for n, crispr in enumerate(crisprs, 1):
      if crispr['Evidence_Level'] >= 3:
         print('>{}_CRISPR_{}_DR'.format(id, n))
         print(crispr['DR_Consensus'])
         left_flank, right_flank = crispr['Regions'][0], crispr['Regions'][-1]
         assert(left_flank['Type'] == 'LeftFLANK')
         assert(right_flank['Type'] == 'RightFLANK')
         print('>{}_CRISPR_{}_leftflank_leader_{}'\
               .format(id, n, left_flank['Leader']))
         print(left_flank['Sequence'])
         print('>{}_CRISPR_{}_rightflank_leader_{}'\
               .format(id, n, right_flank['Leader']))
         print(right_flank['Sequence'])
