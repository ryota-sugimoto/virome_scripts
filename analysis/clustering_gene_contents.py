#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('contig_protein', type=argparse.FileType('r'))
parser.add_argument('contig_labels', type=argparse.FileType('r'))
parser.add_argument('protein_labels', type=argparse.FileType('r'))
args = parser.parse_args()

protein_d = {}
protein_set = []
for s in args.contig_protein:
  l = s.strip().split()
  contig = l[0]
  protein = l[1]
  protein_set.append(protein)
  if contig in protein_d:
    protein_d[contig].append(protein)
  else:
    protein_d[contig] = [protein]

for contig in protein_d:
  protein_d[contig] = list(set(protein_d[contig]))



protein_set = set(protein_set)

import numpy as np
from scipy.sparse import csr_matrix, find
n_rows = len(protein_d.keys())
n_cols = len(protein_set)
contigs = sorted(protein_d.keys())
contigs_index = { contig:i for i,contig in enumerate(contigs) }
proteins = sorted(list(protein_set))
proteins_index = { protein:i for i,protein in enumerate(proteins) }
row = []
col = []
for contig in contigs:
  for protein in protein_d[contig]:
    row.append(contigs_index[contig])
    col.append(proteins_index[protein]) 
data = np.ones(len(row))
init_mat = csr_matrix((data, (row, col)),
                       shape=(n_rows, n_cols), dtype=np.int8)

import pandas
protein_df = pandas.DataFrame(init_mat.todense(),
                                       index=contigs, columns=proteins)



label_l = []
for s in args.contig_labels:
  if s[0] == '#':
    continue
  l = s.strip().split('\t')
  contig,length,gc,db,host,capsid_coding = l[0], int(l[1]), float(l[2]), \
                        list(map(int, l[3:7])), l[16], int(l[18])
  label_l.append((contig,length,gc) + tuple(db) + (host, capsid_coding))
 

label_l.sort()
contigs = [ t[0] for t in label_l ]
contig_label_df =  pandas.DataFrame([t[1:] for t in label_l],
                             index=contigs,
        columns=('length', 'GC', 'NCBI virus', 'NCBI plasmid',
                 'IMG/VR', 'GVD', 'host', 'capsid coding'))

protein_category_cmap = {'DNA_polymerase_family_A': "#56C1FF" ,
                    'DNA_polymerase_family_B': "#01A2FF",
                    'DNA_polymerase_III_subunits': "#0076BA",
                    'DNA_polymerase_family_Y': "#004D80",
                    'RNA_polymerase_subunits': "#16E7CF",
                    'reverse_transcriptase': "#00A89D",
                    'capsid': "#FFFC66",
                    'tail': "#FAE232",
                    'terminase_large_subunit': "#F8BA00",
                    'terminase_small_subunit': "#FF9300",
                    'portal': "#AA7942",
                    'integrase': "#FF644E",
                    'conjugation': "#EF5FA7",
                    'antirestriction': "#CB297B",
                    'NA': "#ffffff" }
label_l = []
for s in args.protein_labels:
 l = s.strip().split()
 protein, label = l[0], l[1]
 if protein in protein_set:
   label_l.append((protein, label))
label_l.sort()
protein_label_df = pandas.DataFrame([[protein_category_cmap[t[1]]
                                       for t in label_l]],
                                    index=['protein_label'],
                                    columns=[ t[0] for t in label_l])
protein_label_df.loc['protein_label'] = protein_label_df.loc['protein_label'].fillna("#ffffff")

df = pandas.concat([protein_df, contig_label_df], axis=1, join='inner')
#df.index.name = 'Contig'
#df_protein_color = pandas.concat([protein_label_df, df], sort=True)

#df_protein_color.loc['protein_label'] = df_protein_color.loc['protein_label'].fillna("#ffffff")

df.loc[df['length'] <= 20000, 'length category'] = 'small'
df.loc[(df['length'] > 20000) & (df['length'] < 200000),
       'length category'] = 'middle'
df.loc[df['length'] >= 200000, 'length category'] = 'large'
length_category_cmap = {'small': '#48C9B0',
                        'middle': '#F5B041',
                        'large': '#C0392B'}
df['length category colors'] = df['length category'].map(length_category_cmap)

df.loc[df['GC'] <= 0.55, 'GC category'] = 'GC low'
df.loc[df['GC'] > 0.55, 'GC category'] = 'GC high'
GC_category_cmap = {'GC low': '#76B4E2',
                    'GC high':'#DD6E43'}
df['GC category colors'] = df['GC category'].map(GC_category_cmap)

db_category_cmap = {0:'#ffffff', 1:'#01A2FF'}
df['ncbi virus colors'] = df['NCBI virus'].map(db_category_cmap)

db_category_cmap = {0:'#ffffff', 1:'#60D836'}
df['ncbi plasmid colors'] = df['NCBI plasmid'].map(db_category_cmap)

db_category_cmap = {0:'#ffffff', 1:'#F9BA00'}
df['imgvr colors'] = df['IMG/VR'].map(db_category_cmap)

db_category_cmap = {0:'#ffffff', 1:'#FF2500'}
df['gvd colors'] = df['GVD'].map(db_category_cmap)

from matplotlib import cm
from matplotlib import colors

'''
jet = cm.get_cmap('gist_rainbow', len(df['host'].unique()))
host_category_cmap = {}
for i, rgb in zip(sorted(df['host'].unique()),
                  jet(range(len(df['host'].unique())))):
  host_category_cmap[i] = colors.rgb2hex(rgb)
host_category_cmap['NA'] = '#ffffff'
'''
host_category_cmap = {'Firmicutes': '#1F77B4',
                      'Bacteroidetes': '#FF7F0E',
                      'Actinobacteria': '#2CA02C',
                      'Proteobacteria': '#D62728',
                      'Fusobacteria': '#9467BD',
                      'Verrucomicrobia': '#8C564B',
                      'Euryarchaeota': '#E377C2',
                      'multiple': '#96989A',
                      'NA': '#ffffff'}
df['host colors'] = df['host'].map(host_category_cmap)

db_category_cmap = {0:'#ffffff', 1:'#FFE24E'}
df['capsid coding colors'] = df['capsid coding'].map(db_category_cmap)

import seaborn as sns
import os

#large
labels = ['ncbi virus colors',
          'ncbi plasmid colors',
          'imgvr colors',
          'gvd colors',
          'length category colors',
          'GC category colors',
          'host colors',
          'capsid coding colors']

sub_df = df.loc[df['length']>=20000]
print(sub_df)
protein_mat = sub_df.loc[:,proteins].astype(int)
top1000_proteins = protein_mat.sum().sort_values().tail(1000).index
top1000_submat = protein_mat.loc[:, top1000_proteins]
labels_subdf = sub_df.loc[:, labels]
large_df = labels_subdf.join(top1000_submat, how='inner')
for lbl in labels:
  protein_label_df.loc[:,lbl] = 'NA'
cdf = pandas.concat([protein_label_df, large_df], join='inner')
cdf.index.name = 'Contig'
with open('large.tsv', 'w') as f:
  f.write(cdf.to_csv(sep='\t'))

sub_df = df.loc[df['length']<20000]
protein_mat = sub_df.loc[:,proteins].astype(int)
top1000_proteins = protein_mat.sum().sort_values().tail(1000).index
top1000_submat = protein_mat.loc[:, top1000_proteins]
labels_subdf = sub_df.loc[:, labels]
large_df = labels_subdf.join(top1000_submat, how='inner')
for lbl in labels:
  protein_label_df.loc[:,lbl] = 'NA'
cdf = pandas.concat([protein_label_df, large_df], join='inner')
cdf.index.name = 'Contig'
with open('small.tsv', 'w') as f:
  f.write(cdf.to_csv(sep='\t'))

sub_df = df.loc[df['length']>0]
protein_mat = sub_df.loc[:,proteins].astype(int)
top1000_proteins = protein_mat.sum().sort_values().tail(1000).index
top1000_submat = protein_mat.loc[:, top1000_proteins]
labels_subdf = sub_df.loc[:, labels]
large_df = labels_subdf.join(top1000_submat, how='inner')
for lbl in labels:
  protein_label_df.loc[:,lbl] = 'NA'
cdf = pandas.concat([protein_label_df, large_df], join='inner')
cdf.index.name = 'Contig'
with open('all.tsv', 'w') as f:
  f.write(cdf.to_csv(sep='\t'))

sub_df = df.loc[df['length']>200000]
protein_mat = sub_df.loc[:,proteins].astype(int)
top1000_proteins = protein_mat.sum().sort_values().tail(2000).index
top1000_submat = protein_mat.loc[:, top1000_proteins]
labels_subdf = sub_df.loc[:, labels]
large_df = labels_subdf.join(top1000_submat, how='inner')
for lbl in labels:
  protein_label_df.loc[:,lbl] = 'NA'
cdf = pandas.concat([protein_label_df, large_df], join='inner')
cdf.index.name = 'Contig'
with open('verylarge.tsv', 'w') as f:
  f.write(cdf.to_csv(sep='\t'))



'''
g = sns.clustermap(top1000_submat,
                   row_colors=[sub_df['ncbi virus colors'],
                               sub_df['ncbi plasmid colors'],
                               sub_df['imgvr colors'],
                               sub_df['gvd colors'],
                               sub_df['GC category colors'],
                               sub_df['host colors'],
                               sub_df['Gram-positive colors']],
                   col_colors=[df_protein_color.loc['protein_label', top1000_proteins]],
                   cmap=['#ffffff', '#834190'])
#                   edgecolor="none")
g.cax.set_visible(False)
g.ax_heatmap.set_xticklabels([])
g.ax_heatmap.set_xticks([])
g.ax_heatmap.set_yticklabels([])
g.ax_heatmap.set_yticks([])
g.savefig('large_contigs.cluster.pdf', dpi=2000)
exit()
#small
sub_df = df.loc[df['length']<20000]
protein_mat = sub_df.loc[:,proteins].astype(int)
top1000_proteins = protein_mat.sum().sort_values().tail(1000).index
top1000_submat = protein_mat.loc[:, top1000_proteins]
g = sns.clustermap(top1000_submat,
                   row_colors=[sub_df['ncbi virus colors'],
                               sub_df['ncbi plasmid colors'],
                               sub_df['imgvr colors'],
                               sub_df['gvd colors'],
                               sub_df['GC category colors'],
                               sub_df['host colors'],
                               sub_df['Gram-positive colors']],
                   col_colors=[df_protein_color.loc['protein_label', top1000_proteins]],
                   cmap=['#ffffff', '#834190'])
#                   edgecolor="none",
g.cax.set_visible(False)
g.ax_heatmap.set_xticklabels([])
g.ax_heatmap.set_xticks([])
g.ax_heatmap.set_yticklabels([])
g.ax_heatmap.set_yticks([])
g.savefig('small_contigs.cluster.pdf', dpi=1000)

#verylarge
sub_df = df.loc[df['length']>=200000]
protein_mat = sub_df.loc[:,proteins].astype(int)
top1000_proteins = protein_mat.sum().sort_values().tail(1000).index
verylarge_submat = protein_mat.loc[:, top1000_proteins]
g = sns.clustermap(verylarge_submat,
                   row_colors=[sub_df['ncbi virus colors'],
                               sub_df['ncbi plasmid colors'],
                               sub_df['imgvr colors'],
                               sub_df['gvd colors'],
                               sub_df['GC category colors'],
                               sub_df['host colors'],
                               sub_df['Gram-positive colors']],
                   col_colors=[df_protein_color.loc['protein_label', top1000_proteins]],
                   cmap=['#ffffff', '#834190'])
#                   edgecolor="none",
#                   alpha=0.7,
#                   linewidths=0)
g.cax.set_visible(False)
g.ax_heatmap.set_xticklabels([])
g.ax_heatmap.set_xticks([])
#g.ax_heatmap.set_yticklabels([])
#g.ax_heatmap.set_yticks([])
g.savefig('verylarge_contigs.cluster.pdf', dpi=1000)

#all
sub_df = df.loc[df['length']>0]
protein_mat = sub_df.loc[:,proteins].astype(int)
top1000_proteins = protein_mat.sum().sort_values().tail(1000).index
top1000_submat = protein_mat.loc[:, top1000_proteins]
g = sns.clustermap(top1000_submat,
                   row_colors=[sub_df['ncbi virus colors'],
                               sub_df['ncbi plasmid colors'],
                               sub_df['imgvr colors'],
                               sub_df['gvd colors'],
                               sub_df['length category colors'],
                               sub_df['GC category colors'],
                               sub_df['host colors'],
                               sub_df['Gram-positive colors']],
                   col_colors=[df_protein_color.loc['protein_label', top1000_proteins]],
                   cmap=['#ffffff', '#834190'])
#                   edgecolor="none",
#                   alpha=0.7,
#                   linewidths=0)
g.cax.set_visible(False)
g.ax_heatmap.set_xticklabels([])
g.ax_heatmap.set_xticks([])
g.ax_heatmap.set_yticklabels([])
g.ax_heatmap.set_yticks([])
g.savefig('all_contigs.cluster.pdf', dpi=1000)
'''
