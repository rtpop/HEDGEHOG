# HEDGEHOG <img src="logo/HEDGEHOG.png" style="height: 1em; vertical-align: middle;">

## Introduction
Community detection in large scale biological networks, such as gene regulatory networks (GRNs) is challenging.

Here we present HEDGEHOG, a tool for community detection in genome wide [PANDA](https://netzoo.github.io/zooanimals/panda/panda/) [1] networks. HEDGEHOG contains two modules: ``bihidef``, an adaptation of [HiDef](https://github.com/fanzheng10/HiDeF) [2] for bi-partite networks, and ``hedgehog`` which provides preocessing and filtering of PANDA GRNs for community detection with BiHidef and processing of BiHidef outputs into "regulatory pathways". 

BiHiDef can be used for any bi-partite network as long as it is formatted as a comma separated edgelist, however HEDGEHOG provides filtering and network processing specific to PANDA networks.

## Installation guide
HEDGEHOG is available as a GitHub python package and can be installed with ``pip``.

```
git clone https://github.com/rtpop/HEDGEHOG.git
cd HEDGEHOG
pip .
```

## Method description

## References
1. Glass K, Huttenhower C, Quackenbush J, Yuan GC. Passing messages between biological networks to refine predicted interactions. PLoS One. 2013 May 31;8(5):e64832. doi: 10.1371/journal.pone.0064832. PMID: 23741402; PMCID: PMC3669401.
2. Zheng, F., Zhang, S., Churas, C. et al. HiDeF: identifying persistent structures in multiscale ‘omics data. Genome Biol 22, 21 (2021). https://doi.org/10.1186/s13059-020-02228-4
