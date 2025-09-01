# HEDGEHOG <img src="logo/HEDGEHOG.png" style="height: 1em; vertical-align: middle;">

## Introduction
Community detection in large scale biological networks, such as gene regulatory networks (GRNs) is challenging.

Here we present HEDGEHOG, a tool for community detection in genome wide [PANDA](https://netzoo.github.io/zooanimals/panda/panda/)<sup>1</sup> networks. HEDGEHOG contains two modules: ``bihidef``, an adaptation of [HiDef](https://github.com/fanzheng10/HiDeF)<sup>2</sup> for bi-partite networks, and ``hedgehog`` which provides preocessing and filtering of PANDA GRNs for community detection with BiHidef and processing of BiHidef outputs into "regulatory pathways". 

BiHiDef can be used for any bi-partite network as long as it is formatted as a comma separated edgelist, however HEDGEHOG provides filtering and network processing specific to PANDA networks.

## Installation guide
HEDGEHOG is available as a GitHub python package and can be installed with ``pip``.

```
git clone https://github.com/rtpop/HEDGEHOG.git
cd HEDGEHOG
pip .
```

## Method description
This tool comprises two modules. The main community detection algorithm is implemented in the ``bihidef`` module. This is an adaptation of Hidef<sup>2</sup> to bi-partite networks. Similar to HiDef, BiHidef identifies network communities that persist across multiple resolutions. To do so, it constructs a proximity graph across multiple resolution values and detects communities at each resolution using [CONDOR](https://netzoo.github.io/zooanimals/comm/condor/)<sup>3</sup> for each node set of the bi-partite network. Communities from adjacent resolutions are considered persistent if their Jaccard similarity exceeds a user-defined threshold. The persistent communities are then organized hierarchically using containment indices, producing a structure of hierarchical, overlapping communities by assigning each lower-level community to its best-matching higher-level counterpart. This is done independently for each of the network partitions, resulting in two higherarchical community structures, one for each node set. A co-clustering relationship is calculated at each resolution based on a similarity score $S_{ij}$ between cluster _i_ in partition 1 and cluster _j_ in partition 2:

```math
S_{ij} = \sum_{u \in C^{(1)}_i} \sum_{v \in C^{(2)}_j} \mathrm{A}_{u,v}
```
<p>$C^{(1)}_{i}$ is the set of nodes in cluster $i$ of partition 1, $C^{(2)}_{j}$ is the set of nodes in cluster $j$ of partition 2, and $A_{u,v}$ is the entry in the adjacency matrix corresponding to the edge (or its weight) between nodes $u$ and $v$.

Due to the nature of bi-partite communities, the co-clustering relationship is not necessarily reciprocal. Meaning that given community $i$ in partition 1 and community $j$ in partition 2, $S_{ij}$ is not necessarily the same as $S_{ji}$. Thus, the co-clustering relationship is calculated for each partition.</p>

The ``hedgehog`` module filters and formates PANDA gene regulatory networks so they are suitable for use with BiHiDef. Furhermore, it also formats and outputs the results of BiHiDef as a ``gmt`` file which can be used as "regulatory pathways" for Gene Set Enrichment Analysis or other downstream applications. 

## Usage

BihiDef can be run on any bi-partite network and takes as input a comma separated edge list. However, highly dense networks are unlikely to yield any workable communities and should be filtered prior to applying BiHidef. BiHiDef will output several files:

* ``pvr.nodes``: communities on the "regulator" partition. This refers to the first partition.
* ``pvg.nodes``: communities on the "target" partition. This refers to the second partition.
* ``pvr.edges`` & ``pvg.edges``: the hierarchical community structure for the two partitions.
* ``cocluster_pvr_Reg.txt`` & ``cocluster_pvg_Tar.txt``: co-clustering relationships for the two nodesets.

### BiHiDef only

```
# import module
from HEDGEHOG import bihidef

# run bihidef
bihidef.bihidef(filename = "edgelist.csv")
```

### BiHiDef & hedgehog on PNADA networks

```
# import library 
from HEDGEHOG import hedgehog
from HEDGEHOG import bihidef

# process panda network
hedgehog.process_edge_list(input_file = "panda_net.txt", output_file = "panda_edgelist.csv")

# run bihidef
bihidef.bihidef(filename = "panda_edgelist.csv")

# Select communities in the gene space based on size
communities = hedgehog.select_communities(filename= "genes.nodes", args.min_size, args.max_size, args.logs)

# output results as gmt
hedgehog.gmt_from_bihidef(
```

## References
1. Glass, K., Huttenhower, C., Quackenbush, J., & Yuan, G. C. (2013). Passing messages between biological networks to refine predicted interactions. PloS one, 8(5), e64832. https://doi.org/10.1371/journal.pone.0064832
2. Zheng, F., Zhang, S., Churas, C. et al. HiDeF: identifying persistent structures in multiscale ‘omics data. Genome Biol 22, 21 (2021). https://doi.org/10.1186/s13059-020-02228-4
3. Platig, J., Castaldi, P. J., DeMeo, D., & Quackenbush, J. (2016). Bipartite Community Structure of eQTLs. PLoS computational biology, 12(9), e1005033. https://doi.org/10.1371/journal.pcbi.1005033
