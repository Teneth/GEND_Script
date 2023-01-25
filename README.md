# GEND_Script
A R language script for calculating gene-gene correlations from single cell RNA sequencing data and then collapsing them into associated networks

Dependencies are the R libraries Rtools, propagate, ggplot2, dplyr

RTools and propagate can be replaced by the basic Cor function, but only if your dataset is small enough to fit your memory limits

GEND- Gene Expression Network Discovery
Attempts to take the upper half of correlated gene space and distribute into patterns of gene expression, networks of similar gene expression. Distinct from GRNs, which attempt to create transcription factor hubs for their networks via various methods. GEND is unbiased and will detect changes in gene expression that are shared across clusters and samples, a facet that standard methods lack.


