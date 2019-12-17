#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('protospacers_bedfile',
                    type=argparse.FileType('r'))
parser.add_argument('-c', '--min_cluster_size',
                    type=int,
                    default=3)
parser.add_argument('-t', '--min_co_occurance',
                    type=int,
                    default=3)
parser.add_argument('-d', '--min_dice_coefficient',
                    type=float,
                    default=0.05)
parser.add_argument('-a', '--community_detection_algorithm',
                    type=str,
                    choices=['labelpropagation', 'walktrap'],
                    default='labelpropagation')
parser.add_argument('-k', '--centrality_algorithm',
                    type=str,
                    choices=['pagerank', 'eigenvector', 'degree'],
                    default='pagerank')
parser.add_argument('-s', '--min_keyspacer_score',
                    type=float,
                    default=0.7)
args = parser.parse_args()

d = {}
for s in args.protospacers_bedfile:
  l = s.strip().split()
  contig = l[0]
  spacer = '|'.join(sorted(l[3].split('|')))
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

import igraph

i,j,v = find(spacer_cooccur_mat)
diag = spacer_cooccur_mat.diagonal()
w = np.where(np.logical_and(2*v/(diag[i]+diag[j]) >= args.min_dice_coefficient,
                            v >= args.min_co_occurance), v, 0)
spacer_cooccur_mat_ = csr_matrix((w, (i,j)),
                                 shape=spacer_cooccur_mat.shape)
spacer_cooccur_mat_.setdiag(0)
spacer_cooccur_mat_.eliminate_zeros()
row, col = spacer_cooccur_mat_.nonzero()
weight = spacer_cooccur_mat_.data
g = igraph.Graph(list(zip(row.tolist(), col.tolist())),
                 vertex_attrs={'name': spacers},
                 edge_attrs={'weight': data.tolist()})
components = g.components()

def cov_d(d, names):
  keys = []
  for name in names:
    keys += list(key for key in d if name in d[key])
  return len(set(keys))

res = []
report = []
for comp in components:
  if len(comp) == 1:
    if args.min_cluster_size == 1:
      res.append(['*' + g.subgraph(comp).vs['name'][0]])
    continue


  comp_g = g.subgraph(comp)
  if args.community_detection_algorithm == 'labelpropagation':
    communities = comp_g.community_label_propagation(weights='weight')
  elif args.community_detection_algorithm == 'walktrap':
    communities = comp_g.community_walktrap(weights='weight').as_clustering()
  else:
    raise('something is wrong')

  for comm in communities:
    N_comm = int(len(comm))
    if not N_comm >= args.min_cluster_size:
      continue
    comm_g = comp_g.subgraph(comm)
    comm_spacer_names = np.array(comm_g.vs['name'])
    report.append((N_comm,
                   comm_g.assortativity_degree(directed=False),
                   comm_g.transitivity_undirected()))
    
    if args.centrality_algorithm == 'eigenvector':
      cscores = np.array(comm_g.evcent(directed=False, weights='weight'))
    elif args.centrality_algorithm == 'pagerank':
      cscores = np.array(comm_g.pagerank(directed=False, weights='weight'))
      cscores /= np.max(cscores)
    elif args.centrality_algorithm == 'degree':
      cscores = np.array(comm_g.degree(loops=False))
      cscores /= np.max(cscores)
    else:
      raise('something is wrong')
    keyspacers_index = (cscores >= args.min_keyspacer_score).nonzero()
    keyspacers = comm_spacer_names[keyspacers_index]

    res.append(list(name if name not in keyspacers else '*'+name for name in comm_spacer_names))
    

for i,((N_comm, assort, clstr), names) in enumerate(zip(report, res)):
  print('\t'.join([str(i), str(N_comm), str(assort), str(clstr), 
                   str(','.join(map(str, names)))]))
