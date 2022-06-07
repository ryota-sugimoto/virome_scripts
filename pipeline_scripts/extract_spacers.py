#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('DR_masked_fastq', type=argparse.FileType('r'))
parser.add_argument('--min_spacer_length', '-s', type=int, default=20)
parser.add_argument('--max_spacer_length', '-l', type=int, default=50)
parser.add_argument('--sample_id', '-n', type=str, default='')
args = parser.parse_args()

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import re
from collections import Counter

re_obj = re.compile('R[^R]+R')

i = 0
dr_containig_seqs = []
for record in SeqIO.parse(args.DR_masked_fastq, 'fastq'):
  for seq in re_obj.findall(str(record.seq)):
    seq = seq[1:-1]
    if len(seq) < args.min_spacer_length or len(seq) > args.max_spacer_length:
      continue
    description = 'fastq_id={}'.format(record.id)
    description += '\t' + '\t'.join(record.description.split('\t')[1:])
    if args.sample_id:
      description += '\tsample_id={}'.format(args.sample_id)
    i += 1
    id = 'SPACER_{}'.format(i)
    spacer_record = SeqRecord(Seq(seq),
                              id=id,
                              name=id,
                              description=description)
    SeqIO.write(spacer_record, sys.stdout, 'fasta')
