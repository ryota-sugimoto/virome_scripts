# Spacer clustering scripts
This directory contains scripts used for the spacer clustering.

## cluster_stats.py
Calculate clustering coefficient from clustering result. The script requires protospacer BED file and clustering result. The output is seperated by tabular and contains three columns; the first column is the cluster ID, the second column is the clustering coefficient, and the third column is the list of spacers in the cluster.

### Example
```
./cluster_stats.py all_protospacers.bed protospacers_cluster.tab > cluster_stats
```

## generate_abc_edgefile.py
Calculate a co-occurrence graph from initially clustered protospacer bed file and output the graph in the abc format. The bed file must be initially clustered based on the distance (see the protocol).

### Example
```
./generate_abc_edgefile.py initial_cluster.bed > protospacers.abc
```

## mark_crispr_targeted.py
Mark CRISPR targeted regions based on the spacer clusters. The script requires the original protospacer BED file and the spacer clustering result (see the protocol). This script calculates the median distances between the protospacers from each spacer cluster, and marks the regions around the protospacers according to this distances.

### Example
```
./mark_crispr_targeted.py all_protospacers.bed protospacers_cluster.tab > crispr_targeted.bed
```
