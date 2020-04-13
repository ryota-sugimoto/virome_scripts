#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('protospacers_bedfile',
                    type=argparse.FileType('r'))
parser.add_argument('-t', '--min_co_occurance',
                    type=int,
                    default=3)
parser.add_argument('-d', '--min_dice_coefficient',
                    type=float,
                    default=0.3)
args = parser.parse_args()

d = {}
for s in args.protospacers_bedfile:
  l = s.strip().split()
  contig = l[0]
  spacer = l[3]
  if contig in d:
    d[contig].add(spacer)
  else:
    d[contig] = set([spacer])

from itertools import chain
import numpy as np
np.set_printoptions(precision=3)
spacers = np.array(sorted(list(set(chain(*d.values())))))
spacers_id = {key:i for i,key in enumerate(spacers)}

contigs = np.array(sorted(d.keys()))
contigs_id = {key:i for i,key in enumerate(contigs)}

row = []
col = []
for contig in contigs:
  for spacer in d[contig]:
    row.append(contigs_id[contig])
    col.append(spacers_id[spacer])

data = np.ones(len(row))

from scipy.sparse import csr_matrix, find
contig_spacer_mat = csr_matrix((data, (row,col)),
                               shape=(len(contigs), len(spacers)))

spacer_cooccur_mat = contig_spacer_mat.T * contig_spacer_mat

i,j,v = find(spacer_cooccur_mat)
diag = spacer_cooccur_mat.diagonal()

w = np.where(np.logical_and(2*v/(diag[i]+diag[j]) >= args.min_dice_coefficient,
                            v >= args.min_co_occurance), v, 0)
spacer_cooccur_mat_ = csr_matrix((w, (i,j)),
                                 shape=spacer_cooccur_mat.shape)
spacer_cooccur_mat_.setdiag(0)
spacer_cooccur_mat_.eliminate_zeros()
from scipy.sparse import triu

for i,j,v in zip(*find(triu(spacer_cooccur_mat_, k=1))):
  print(spacers[i], spacers[j], v)
