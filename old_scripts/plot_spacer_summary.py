#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('sam', type=argparse.FileType('r'))
args = parser.parse_args()

import numpy as np

v = []
for s in args.sam:
  if s[0] != '@':
    l = s.strip().split()
    spacer_contig = l[0].split('_')[0]
    crispr_id = l[0].split('_')[3]
    cas_type = l[0].split('_')[6]
    protospacer_contig = l[2].split('_')[0]
    v.append((spacer_contig, crispr_id, cas_type, protospacer_contig))

crisprs = list(set(map(lambda t: t[0]+'_'+t[1]+'_'+t[2], v)))
v = np.array(v)
cas_types = np.unique(v[:, 2])

count = []
for crispr in crisprs:
  spacer_contig, crispr_id, cas_type = crispr.split('_')
  w = v[np.logical_and(v[:,0]==spacer_contig,
                       v[:,1]==crispr_id,
                       v[:,2]==cas_type)]
  num_all_spacers = len(w)
  num_aligned_spacers = len(w[w[:, 3]!='*'])
  count.append((crispr, cas_type, num_all_spacers, num_aligned_spacers))
count = np.array(count)


import matplotlib.pyplot as plt
import matplotlib as mpl

color_mapper = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0,
                                                         vmax=len(cas_types)),
                               cmap='jet')

markers = ['+', '^', '*', 'p']
import re
for i,cas_type in enumerate(sorted(cas_types, key=lambda v: v.count('I'))):
  w = count[count[:,1]==cas_type]
  x,y = w[:,2],w[:,3]
  if cas_type == 'None':
    marker = '.'
  else:
    marker = markers[cas_type.count('I')]
  plt.scatter(x, y, c=color_mapper.to_rgba(i),
              s=100, edgecolors='none', alpha=0.6, marker=marker, 
              label=cas_type)
axes = plt.gca()
x_vals = np.array(axes.get_xlim())
plt.plot(x_vals,x_vals, '--')
plt.xlabel('#spacers')
plt.ylabel('#protospacers')
plt.legend(loc='upper left')
plt.show()
