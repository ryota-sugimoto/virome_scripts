#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('protospacers_bedfile',
                    type=argparse.FileType('r'))
parser.add_argument('cluster_file',
                    type=argparse.FileType('r'))
parser.add_argument('-t', '--min_co_occurance',
                    type=int,
                    default=3)
parser.add_argument('-d', '--min_dice_coefficient',
                    type=float,
                    default=0.3)
parser.add_argument('-m', '--min_cluster_size',
                    type=int,
                    default=10)
parser.add_argument('-s', '--min_keyspacer_score',
                    type=float,
                    default=0.7)
args = parser.parse_args()

clusters = []
for s in args.cluster_file:
  clusters.append(s.strip().split())

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

import igraph
from scipy.sparse import triu
upper = triu(spacer_cooccur_mat_)
row, col = upper.nonzero()
weight = upper.data

g = igraph.Graph(list(zip(row.tolist(), col.tolist())),
                 vertex_attrs={'name': spacers},
                 edge_attrs={'weight': weight})

assortativity = g.assortativity_degree(directed=False)

for cluster_id, cluster in enumerate(clusters):
  if len(cluster) < args.min_cluster_size:
    continue
  subg = g.subgraph(cluster)
  clustering_coefficient = subg.transitivity_undirected()
  
  degrees = np.array(subg.strength(loops=False, weights='weight'),
                     dtype=np.float)
  degrees /= np.max(degrees)
  subg_spacers = np.array(subg.vs['name'])
  keyed_subg_spacers = list(s+':'+str(d) 
                       for s,d in zip(subg_spacers, degrees))
  tmp = ','.join(keyed_subg_spacers) 
  print("\t".join(list(map(str,[cluster_id,
                           clustering_coefficient]))) + '\t' + tmp)
