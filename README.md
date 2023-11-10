# Paree_map

This repository contains the scripts and data corresponding to the article XXXX DOI XXXX.

Directories:

trimming: contains the script used to verify the presence of the Bcg1 site in 2b-Rad sequence and trim them to the 36 bp inserted sample DNA fragment

simulations: contains the script used to simulate RIAIL panel construction

maps: contains the estimated genetic position from the RIAILs with the R/qtl package for both rec-1 alleles, including the different intercross types. As well as the code to plot them against the maps from Rockman & Kruglyak, 2009 https://doi.org/10.1371/journal.pgen.1000419

nondisjunction: contains the data of the assay done of non-disjunction (fertility, male proportion, egg viability).

genotype: contains the genotype of the RIAIL used to estimate the linkage maps. Note that the .fastq files of the parental strains are available on ncbi (PRJNA1037511). 

files:

gini&AUC.R: the R functions used to calculate the Lorenz curves, gini coefficients, and the relative Area Under the Curve (AUC)

utils.R: some useful R function


