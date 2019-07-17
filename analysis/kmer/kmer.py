#!/usr/bin/env python3
from collections import defaultdict
from itertools import product
from Bio.Seq import Seq

def reverse_complement(seq):
  comp = {'A': 'T',
          'C': 'G',
          'G': 'C',
          'T': 'A'}
  seq = list(seq)
  seq.reverse()
  return ''.join(map(lambda n:comp[n], seq))

class KmerCounter():
  def __init__(self, k):
    self.k = k
    self._calc_kmer()

  def _calc_kmer(self):
    self.kmer = [Seq(''.join(mer)) for mer in product('ACGT', repeat=self.k)]
    self.kmer.sort()
    pass_set = set([])
    self.balanced_kmer = []
    self.reverse_complement = {}
    for mer in self.kmer:
      comp = mer.reverse_complement()
      self.reverse_complement[mer] = comp
      if mer not in pass_set:
        self.balanced_kmer.append(mer)
      pass_set.add(comp)
      pass_set.add(mer)
    self.balanced_kmer.sort()
 
  def balance_strand(self, count):
    for mer in self.balanced_kmer:
      comp = self.reverse_complement[mer]
      if mer != comp:
        new_value = count[mer] + count[comp]
        count[mer] = count[comp] = new_value
    return count

  def count(self, seq):
    seq = seq.upper()
    count = defaultdict(int)
    for i in range(len(seq)-self.k+1):
      count[seq[i:i+self.k]] += 1
    count = self.balance_strand(count)
    v = []
    for mer in self.balanced_kmer:
      v.append(count[mer])
    return v
       
from Bio import SeqIO
if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument('k', type=int)
  parser.add_argument('fasta', type=argparse.FileType('r'))
  parser.add_argument('--frequency', '-f',
                      action='store_true',
                      default=False)
  args = parser.parse_args()
  
  counter = KmerCounter(args.k)
 
  print('\t'.join(['#seq_id'] + list(map(str,counter.balanced_kmer))))
  for record in SeqIO.parse(args.fasta, 'fasta'):
    count = counter.count(record.seq)
    if args.frequency:
      kmer_count = map(lambda x: x/sum(count), count)
    else:
      kmer_count = count
    print('\t'.join([record.id]+list(map(str,kmer_count))))
