#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('blastfmt6', type=argparse.FileType('r'))
parser.add_argument('contigs_fasta', type=argparse.FileType('r'))
args = parser.parse_args()

from Bio import SeqIO
record_dict = SeqIO.to_dict(SeqIO.parse(args.contigs_fasta, "fasta"))

import sys
d = {}
for s in args.blastfmt6:
  l = s.strip().split()
  contig = l[1]
  
  qstart, qend = int(l[6]), int(l[7])
  sstart, send = int(l[8]), int(l[9])
  
  if contig not in d:
    d[contig] = (qstart, qend, sstart, send)
  else:
    if abs(d[contig][1] - d[contig][0]) < abs(qend - qstart):
      d[contig] = (qstart, qend, sstart, send)
  
for contig in d:
  qstart, qend, sstart, send = d[contig]
  if sstart > send:
    record_dict[contig].seq = record_dict[contig].seq.reverse_complement()
    sstart = len(record_dict[contig].seq) - sstart + 1
    send = len(record_dict[contig].seq) - send + 1
  
  chop_pos = sstart-qstart+1 
  record_dict[contig].seq = record_dict[contig].seq[chop_pos:] \
               + record_dict[contig].seq[:chop_pos] 

  SeqIO.write(record_dict[contig], sys.stdout, "fasta")
