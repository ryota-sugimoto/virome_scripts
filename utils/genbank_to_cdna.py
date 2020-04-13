#!/usr/bin/env python3

import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re

import argparse 

parser = argparse.ArgumentParser()
parser.add_argument('genbank', type=argparse.FileType('r'))
parser.add_argument('-r', '--only_ribosomal_proteins', default=False,
                    action='store_true')
args = parser.parse_args()

gb = SeqIO.parse(args.genbank, 'genbank')

expr = re.compile('(50|30)S ribosomal( subunit)? protein')


for record in gb:
  i = 0
  for feature in record.features:
    if not feature.type == "CDS":
      continue
    if 'pseudo' in feature.qualifiers:
      continue
    if args.only_ribosomal_proteins and \
       not any(map(lambda s: expr.match(s), feature.qualifiers['product'])):
      continue
    if 'parts' not in feature.location.__dict__:
      start = feature.location._start
      end = feature.location._end
      if feature.location._strand == 1:
        seq = record.seq[start:end]
      else:
        seq = record.seq[start:end].reverse_complement()
    else:
      l = []
      for part in feature.location.parts:
        start = part._start
        end = part._end
        strand = part._strand
        if strand == 1:
          l.append(record.seq[start:end])
        else:
          l.append(record.seq[start:end].reverse_complement())
      joined_sequence = ''
      for s in l:
        joined_sequence += s
      seq = joined_sequence
    r = SeqRecord(seq,
                  id=record.id + "_" + str(i),
                  description=feature.qualifiers['product'][0])
    i += 1
    sys.stdout.write(r.format('fasta'))
