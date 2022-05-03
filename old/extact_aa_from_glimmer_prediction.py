#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
parser.add_argument('glimmer_predict', type=argparse.FileType('r'))
parser.add_argument('-n', '--nucleotide', action='store_true', default=False)
args = parser.parse_args()
import sys
from Bio import SeqIO
record_dict = SeqIO.index(args.fasta.name, 'fasta')

from Bio.SeqRecord import SeqRecord

for s in args.glimmer_predict:
  if s[0] == '>':
    seq_id = s.strip().split()[0][1:]
    continue
  
  l = s.strip().split()
  gene_id = l[0]
  gene_start, gene_end, orf, score = int(l[1]), int(l[2]), \
                                     int(l[3]), float(l[4])
  seq_len = len(record_dict[seq_id])
  seq = record_dict[seq_id].seq
  if gene_start > gene_end and orf > 0:
    gene_seq = seq[gene_start-1:] + seq[:gene_end]
  elif gene_start < gene_end and orf < 0:
    gene_seq = seq[gene_end-1:] + seq[:gene_start]
    gene_seq = gene_seq.reverse_complement()
  elif gene_start < gene_end and orf > 0:
    gene_seq = seq[gene_start-1:gene_end]
  elif gene_start > gene_end and orf < 0:
    gene_seq = seq[gene_end-1:gene_start]
    gene_seq = gene_seq.reverse_complement()
  else:
    raise('something is wrong')

  out_id = '{}_{}_{}_{}'.format(seq_id, gene_start-1, gene_end, orf)
  if not args.nucleotide:
    record = SeqRecord(gene_seq.translate(), id=out_id, description='')
  else:
    record = SeqRecord(gene_seq, id=out_id, description='')
    
  SeqIO.write(record, sys.stdout, 'fasta')
