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
    Cas = sequence['Cas']
    loci = list( (cas['Type'], cas['Start']) for cas in Cas)
    crisprs = sequence['Crisprs']
    for n, crispr in enumerate(crisprs, 1):
      name = crispr['Name']
      if crispr['Evidence_Level'] >= 3:
         start = crispr['Start']
         if loci:
           closest_loci = sorted(loci,
                                key=lambda locus: abs(start-locus[1]))[0]
           closest_cas, distance = closest_loci[0], abs(start-closest_loci[1])
         else:
           closest_cas, distance = 'None', 'NA'
         print('>{}#CRISPR_DR#name:{}#Cas:{}#distance:{}'.format(id, name, 
                                                     closest_cas, distance))
         print(crispr['DR_Consensus'])
         left_flank, right_flank = crispr['Regions'][0], crispr['Regions'][-1]
         assert(left_flank['Type'] == 'LeftFLANK')
         assert(right_flank['Type'] == 'RightFLANK')
         print('>{}#CRISPR_leftflank#name:{}#leader:{}'\
               .format(id, n, left_flank['Leader']))
         print(left_flank['Sequence'])
         print('>{}#CRISPR_rightflank#name:{}#leader:{}'\
               .format(id, n, right_flank['Leader']))
         print(right_flank['Sequence'])
