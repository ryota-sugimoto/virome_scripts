#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO
import os
prefix = os.path.splitext(args.fasta.name)[0]

record_dict = SeqIO.index(args.fasta.name, 'fasta')
ids = list(record_dict.keys())


for n,id in enumerate(ids, 1):
  with open("{}_{}_.fasta".format(prefix, n), 'w') as remain_f:
    for id_ in ids:
      print(id,id_)
      if id_ == id:
        with open("{}_{}.fasta".format(prefix, n), 'w') as selected_f:
          SeqIO.write(record_dict[id_], selected_f, 'fasta')
      else:
          SeqIO.write(record_dict[id_], remain_f, 'fasta')
