#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('tetramer', type=argparse.FileType('r'))
args = parser.parse_args()

import numpy as np


data = []
for s in args.tetramer:
  if s[0] == '#': continue
  l = s.strip().split('\t')
  data.append((l[0],l[1],l[2:],[0,0]))
data = np.array(data, dtype=[('name', 'S50'),
                             ('category', 'S50'),
                             ('tetramer', 'f4', (136,)),
                             ('tsne_embed', 'f4', (2,))])
categories = sorted(list(set(data['category'])))
category_y = { cat: y for y,cat in enumerate(categories) }

from sklearn.preprocessing import normalize
data['tetramer'] = normalize(data['tetramer'], norm='l1')

#from sklearn.manifold import TSNE
from openTSNE.sklearn import TSNE
from matplotlib import cm
import matplotlib.pyplot as plt
cmap = cm.get_cmap('jet', len(categories))

for i in [1]:
  perp = i*1000
  tsne = TSNE(perplexity=perp,
              n_jobs=10,
              initialization="pca",
              metric="cosine",
              n_iter=30000,
              neighbors="approx",
              negative_gradient_method="fft")
  data['tsne_embed'] = tsne.fit_transform(data['tetramer'])
  
  for category in categories:
    if str(category.decode('utf-8')) == 'putative_phage':
      marker = 'x'
      size = 0.8
    else:
      marker = None
      size = 0.5
    plt.scatter(data[data['category']==category]['tsne_embed'][:,0],
                data[data['category']==category]['tsne_embed'][:,1],
                marker = marker,
                s=size,
                c=[cmap(category_y[category])],
                alpha=0.4,
                label=category.decode('utf-8'))
  plt.legend()
  plt.show()
  plt.savefig('test_perplexity_{}_iter_30000_pca.png'.format(perp),
              dpi=1000)
  plt.close()
