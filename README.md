# HEDGEHOG <img src="logo/HEDGEHOG.png" style="height: 1em; vertical-align: middle;">

## Introduction
Community detection in large scale biological networks, such as gene regulatory networks (GRNs) is challenging.

Here we present HEDGEHOG, a tool for community detection in genome wide PANDA (cite) networks. HEDGEHOG contains two modules: bihidef, an adaptation of HiDef (cite) for bi-partite networks, and hedgehog which provides preocessing and filtering of PANDA GRNs for community detection with BiHidef and processing of BiHidef outputs into "regulatory pathways". 

BiHiDef can be used for any bi-partite network as long as it is formatted as a comma separated edgelist, however HEDGEHOG provides filtering and network processing specific to PANDA networks.

## Installation guide
HEDGEHOG is available as a GitHub python package and can be installed with ``pip``.

```
git clone https://github.com/rtpop/HEDGEHOG.git
cd HEDGEHOG
pip .
```

## Method description
