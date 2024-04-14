# iimmsci2s

<img src="doc/fig-2s-trees.pdf" width="48">

A C program to calculate asymptotic maximum likelihood estimates of parameters under the MSCI model (a) of two species given that the true model is MSCM (b-d), or vice versa, by minimizing the KL divergence.  Assume one sequence per species, an infinite number of loci, and a Poisson model of mutations.  Sequence length can be finite (n) or infinite.

Adapted from a program IMMSci2s written by [Ziheng Yang](http://abacus.gene.ucl.ac.uk/) for [Jiao et al, 2020](https://doi.org/10.1093/sysbio/syaa001).

References:

* [Huang et al. (2022)](https://doi.org/10.1093/molbev/msac237) Inference of gene flow between species under misspecified models. Mol. Biol. Evol., 39(12):msac237
