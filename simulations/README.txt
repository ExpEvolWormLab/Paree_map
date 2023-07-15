# Simulations were used to test the effect of selection and drift on allele frequency deviation
# and to obtain a null distribution of uniform maps with uniform recombination (flat landscape)

domain_boundaries.Rdata: contains the domain limit from Rockman and Kruglyak, 2009 https://doi.org/10.1371/journal.pgen.1000419
=> used to get the chromosome size and extrapolate

simulations_functions.R: script contains general function to set up a Wright-Fisher simulations

riails_simu.R:  script to simulate RIAIL subpanels construction with modified simulation function

simu_map_construction.R: script to estimate genetic linkage map from the simulated RIAILs subpanels

parameters: directory contains .Rdata files with parameters for simulations: 

- num.loci: the number of loci to simulate
- recombination.rate: estimated from consensus genetic linkage map of rec-1
- male.allele: identifier for the the O sexual chromosome for sex determination
- nextinction: the number of line lost during derivation
- ngen.intercross and ngen.selfing the number of generation of intercross and selfing to simulate
- nplates: the number of metapopulations (=petri dishes were the crosses were done)
- sex.chromosome: which loci are the X chromosome (for sex determination and is treated differently because hemyzygosity)

=> one parameters file per cross subpanels (from A to H; details in Table S2 of the associated article)
