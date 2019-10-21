#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('protospacers_bedfile',
                    type=argparse.FileType('r'))
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

min_sig = 0.05
i,j,v = find(spacer_cooccur_mat)
diag = spacer_cooccur_mat.diagonal()
w = np.where( 2*v/(diag[i]+diag[j]) > min_sig, v, 0)
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

print('all_assort', g.assortativity_degree(directed=False))
print('all_transitivity', g.transitivity_undirected())
def cov_d(d, names):
  keys = []
  for name in names:
    keys += list(key for key in d if name in d[key])
  return len(set(keys))

N = len(d)
print('N', N)
res = []
for comp in components:
  if len(comp) == 1:
    res.append(['*' + g.subgraph(comp).vs['name'][0]])
    continue

  comp_g = g.subgraph(comp)
  print('N_comp', len(comp))
  print('comp_assort', comp_g.assortativity_degree(directed=False))
  print('comp_transitivity', comp_g.transitivity_undirected())
  #communities = comp_g.community_label_propagation(weights='weight')
  #communities = comp_g.community_multilevel(weights='weight')
  #communities = comp_g.community_spinglass(weights='weight')
  communities = comp_g.community_walktrap(weights='weight').as_clustering()
  
  for comm in communities:
    N_comm = len(comm)
    comm_g = comp_g.subgraph(comm)
    comm_spacer_names = np.array(comm_g.vs['name'])
    print('N_comm', N_comm)
    print('comm_assort', comm_g.assortativity_degree(directed=False))
    print('comm_transitivity', comm_g.transitivity_undirected())
    
    '''
    cscores = np.array(comm_g.evcent(directed=False, weights='weight'))
    print('eigenvalue', cscores, np.sum(cscores))
    keyspacers_index = (cscores >= 0.9).nonzero()
    keyspacers = set(comm_spacer_names[keyspacers_index])
    print('eigenvalue', N_comm, len(keyspacers_index[0]), cov_d(d, keyspacers))
    '''
 
    cscores = np.array(comm_g.pagerank(directed=False, weights='weight'))
    #print('pagerank', cscores/np.max(cscores), np.sum(cscores))
    cscores /= np.max(cscores)
    keyspacers_index = (cscores >= 0.9).nonzero()
    keyspacers = comm_spacer_names[keyspacers_index]
    #print('pagerank', N_comm, len(keyspacers_index[0]), cov_d(d, keyspacers))

    '''
    cscores = np.array(comm_g.degree(loops=False))
    #keyspacers_index = np.argsort(cscores)[-min(N_comm, top_n):]
    keyspacers_index = (cscores >= 0.8 * np.max(cscores)).nonzero()
    keyspacers = comm_spacer_names[keyspacers_index]
    print('degree', N_comm, len(keyspacers_index[0]), cov_d(d, keyspacers))
    '''

    res.append(list(name if name not in keyspacers else '*'+name for name in comm_spacer_names))

for names in res:
  print(','.join(map(str,names)))
