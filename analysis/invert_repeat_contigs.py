#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys

parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=argparse.FileType('r'))
parser.add_argument('--seed_length', type=int, default=23)
parser.add_argument('--slide_size', type=int, default=10)
parser.add_argument('--num_slide', type=int, default=30)
parser.add_argument('--min_iterate_alignment_score', type=float, default=0.9)
parser.add_argument('--max_iterate', type=int, default=100)
parser.add_argument('--max_fail_count', type=int, default=10)
parser.add_argument('--max_alignment_length', type=int, default=1000)
parser.add_argument('--min_alignment_length', type=int, default=40)
parser.add_argument('--show_alignment', default=False,
                    action='store_true')
parser.add_argument('--no_trim' , default=False,
                    action='store_true')
args = parser.parse_args()

for record in SeqIO.parse(args.fasta, 'fasta'):
  for i in range(args.num_slide):
    seed_begin = i*args.slide_size
    seed_end = seed_begin + args.seed_length
    seed = record.seq[seed_begin:seed_end]
    if len(seed) != args.seed_length:
      break
    seed_rc = seed.reverse_complement()
    pos = record.seq.rfind(seed_rc,
                           int(len(record.seq)*0.9))
    if pos >= 0 and 'N' not in seed:
      offset = 0
      expand = 1
      fail_count = 0
      last_alignment_length = -1
      for _ in range(args.max_iterate):
        alignment_length = args.seed_length + offset + expand
        if fail_count > args.max_fail_count \
           or alignment_length > args.max_alignment_length \
           or alignment_length == last_alignment_length:
          break
        seq1 = record.seq[seed_begin:seed_begin+alignment_length]
        seq2 = record.seq[pos+args.seed_length-alignment_length:pos+args.seed_length]
        seq2_rc = seq2.reverse_complement()
        if args.show_alignment:
          alignments = pairwise2.align.globalms(seq1, seq2_rc, 1, -1, -1, -1)
          score = alignments[0][2]
        else:
          score = pairwise2.align.globalms(seq1, seq2_rc, 1, -1, -1, -1,
                                           score_only=True)
        if score > alignment_length * args.min_iterate_alignment_score:
          offset = offset + expand
          expand *= 2
        else:
          fail_count += 1
          expand = 1
        last_alignment_length = alignment_length
      if args.show_alignment:
        print('ir_begin',seed_begin)
        print('seq_len', len(record.seq))
        print('seed_pos', pos+args.seed_length)
        print('alignment_length', alignment_length)
        print(format_alignment(*alignments[0]))
      if score > alignment_length * 0.8 \
             and alignment_length > args.min_alignment_length:
        if not args.no_trim:
          record.seq = record.seq[seed_begin:pos+args.seed_length] 
        SeqIO.write(record, sys.stdout, "fasta")
        sys.stdout.flush()
      break
