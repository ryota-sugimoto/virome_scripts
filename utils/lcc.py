#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO
from Bio.SeqUtils import lcc

for record in SeqIO.parse(args.fasta, 'fasta'):
  complexity = lcc.lcc_simp(record.seq)
  print("{}\t{}".format(record.id, complexity))
