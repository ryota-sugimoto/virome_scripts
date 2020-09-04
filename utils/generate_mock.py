#!/usr/bin/env python3

import sys
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO
from collections import Counter
from random import choices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

nuc = 'ATGC'

for record in SeqIO.parse(args.fasta, 'fasta'):
  c = Counter(record.seq)
  id = record.id + '_mock'
  seq = ''.join(choices(nuc, weights=[c[n] for n in nuc], k=len(record.seq)))
  SeqIO.write(SeqRecord(Seq(seq),
                        id=id,
                        description=''),
              sys.stdout,
              'fasta')
