#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('contig_protospacers',
                    type=argparse.FileType('r'))
args = parser.parse_args()

d = {}
for s in args.contig_protospacers:
  l = s.strip().split()
  d[l[0]] = set(l[1].split(','))

from itertools import chain
spacers = set(chain(*d.values()))
spacers_id = {key:i for i,key in enumerate(sorted(list(spacers)))}

contigs = d.keys()
contigs_id = {key:i for i,key in enumerate(sorted(list(contigs)))}

row = []
col = []
for contig in contigs:
  for spacer in d[contig]:
    row.append(contigs_id[contig])
    col.append(spacers_id[spacer])

import numpy as np
np.set_printoptions(precision=3)
data = np.ones(len(row))

from scipy.sparse import csr_matrix, find
from scipy.sparse.csgraph import connected_components
contig_spacer_mat = csr_matrix((data, (row,col)),
                               shape=(len(contigs), len(spacers)))

spacer_cooccur_mat = contig_spacer_mat.T * contig_spacer_mat

#TODO too slow
'''
nj = spacer_cooccur_mat.diagonal()
ni = nj.reshape(-1,1)
dice_coefficients = 2*spacer_cooccur_mat/(ni + nj)
spacer_cooccur_mat[(dice_coefficients<min_sig).nonzero()] = 0
spacer_cooccur_mat.eliminate_zeros()
'''

min_sig = 0.05
i,j,v = find(spacer_cooccur_mat)
diag = spacer_cooccur_mat.diagonal()
w = np.where( 2*v/(diag[i]+diag[j]) > min_sig, v, 0)
spacer_cooccur_mat_ = csr_matrix((w, (i,j)),
                                 shape=spacer_cooccur_mat.shape)
spacer_cooccur_mat_.eliminate_zeros()

from scipy.stats import describe
n_components, labels = connected_components(spacer_cooccur_mat_, directed=False)
print('n_components', n_components)
unique, counts = np.unique(labels, return_counts=True)


N = len(d.keys())
for label,count in zip(unique,counts):
  if count == 1:
    continue
  print('component_label', label)
  print('component_size', count)
  index = (labels==label).nonzero()[0]
  n = spacer_cooccur_mat_[index,:][:,index].todense()
  n_diag = n.diagonal()
  ni = n_diag.reshape(-1,1)
  nj = n_diag
  ll = np.multiply(N, np.log(N))
  ll -= np.multiply(ni, np.ma.log(ni).filled(0)) \
      + np.multiply(nj, np.ma.log(nj).filled(0))
  ll += np.multiply(n, np.ma.log(n).filled(0))
  
  tmp = N - ni - nj + n
  ll += np.multiply(tmp, np.ma.log(tmp).filled(0))
  
  tmp = ni - n
  ll += np.multiply(tmp, np.ma.log(tmp).filled(0))
  
  tmp = nj - n
  ll += np.multiply(tmp, np.ma.log(tmp).filled(0))
  
  tmp1 = N - ni
  tmp2 = N - nj
  ll -= np.multiply(tmp1, np.ma.log(tmp1).filled(0)) \
      + np.multiply(tmp2, np.ma.log(tmp2).filled(0))
  
  ll *= 2
  #ll = csr_matrix(ll)
  print(ll)



'''
import networkx as nx
g = nx.from_scipy_sparse_matrix(spacer_cooccur_mat)
print('size', len(g))
print('n_component', nx.number_connected_components(g))
print('assortability', nx.degree_assortativity_coefficient(g))
for c in nx.connected_components(g):
  subg = g.subgraph(c)
  print('size', len(c))
  print('assortability', nx.degree_assortativity_coefficient(subg))
  print('clustering', nx.average_clustering(subg))
'''

'''
import matplotlib.pyplot as plt
f = plt.figure()
nx.draw(g)
plt.show()
f.savefig('graph.pdf')
plt.show()
'''
