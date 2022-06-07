#!/usr/bin/env python
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('all_protospacers_bed', type=argparse.FileType('r'))
parser.add_argument('protospacers_cluster_tab', type=argparse.FileType('r'))
parser.add_argument('-m', '--min_cluster_size', type=int, default=10)
parser.add_argument('-x', '--max_extend_size', type=int, default=10000)
parser.add_argument('-c', '--min_cluster_ratio', type=float, default=0.3)
args = parser.parse_args()

spacer_clusterid = {}
clusterid_size = {}
for i,s in enumerate(args.protospacers_cluster_tab):
   l = s.strip().split('\t')
   if len(l) < args.min_cluster_size:
     continue
   for spacer in l:
     spacer_clusterid[spacer] = i
   clusterid_size[i] = len(l)

import portion as P
bed_db = {}
begin_positions = {}
for s in args.all_protospacers_bed:
  l = s.strip().split('\t')
  contig = l[0]
  begin = int(l[1])
  end = int(l[2])
  spacer = l[3].split(' ')[0]
  if spacer not in spacer_clusterid:
    continue
  clusterid = spacer_clusterid[spacer]

  if contig in bed_db:
    if clusterid in bed_db[contig]:
      bed_db[contig][clusterid].append(P.closedopen(begin, end))
    else:
      bed_db[contig][clusterid] = [P.closedopen(begin, end)]
  else:
    bed_db[contig] = { clusterid: [P.closedopen(begin, end)] }

  if clusterid in begin_positions:
    if contig in begin_positions[clusterid]:
      begin_positions[clusterid][contig].append(begin)
    else:
      begin_positions[clusterid][contig] = [begin]
  else:
    begin_positions[clusterid] = { contig: [begin] }

begin_diffs = {}
import numpy as np
for clusterid in begin_positions:
  l = []
  for contig in begin_positions[clusterid]:
    l.extend(np.diff(sorted(begin_positions[clusterid][contig])))
  begin_diffs[clusterid] = l

clusterid_extendsize = {}
for clusterid in begin_diffs:
  if len(begin_diffs[clusterid]) < 1:
    clusterid_extendsize[clusterid] = 0
    continue
  extend_size = int(3.45 * np.median(begin_diffs[clusterid])/np.log(2))
  '''
  This is about half the length of 1e-3 right side 
  percentile of the exponential distribution
  '''

  if extend_size > args.max_extend_size:
    clusterid_extendsize[clusterid] = args.max_extend_size
  else:
    clusterid_extendsize[clusterid] = extend_size

merged_bed = {}
for contig in bed_db:
  new_intervals = []
  for clusterid in bed_db[contig]:
    extend = clusterid_extendsize[clusterid]
    intervals = sorted(bed_db[contig][clusterid], key=lambda x: x.lower)
    prev = intervals[0]
    count = 1
    for i in range(1,len(intervals)):
      if prev.upper+extend > intervals[i].lower-extend:
        prev = P.closedopen(prev.lower, intervals[i].upper)
        count += 1
      else:
        if count > clusterid_size[clusterid]*args.min_cluster_ratio:
          if prev.lower - extend < 0:
            lower = 0
          else:
            lower = prev.lower - extend
          upper = prev.upper + extend
          new_interval = P.closedopen(lower, upper)
          new_intervals.append((new_interval, clusterid, count))
        prev = intervals[i]
        count = 1
    if count > clusterid_size[clusterid]*args.min_cluster_ratio:
      if prev.lower - extend < 0:
        lower = 0
      else:
        lower = prev.lower - extend
      upper = prev.upper + extend
      new_interval = P.closedopen(lower, upper)
      new_intervals.append((new_interval, clusterid, count))
  if new_intervals:
    merged_bed[contig] = new_intervals

for contig in sorted(merged_bed.keys()):
  for interval, clusterid, count in merged_bed[contig]:
    print(contig, interval.lower, interval.upper, 
        '{}:{}'.format(clusterid, count), sep='\t')
        
