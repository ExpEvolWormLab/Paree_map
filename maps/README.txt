###
rockman_recombinationmaps.csv: genetic position from reference maps (Rockman & Kruglyak, 2009 https://doi.org/10.1371/journal.pgen.1000419)

###
recmaps: the intercross specific and consensus rec-1 wild-type and mutant genetic linkage maps of this study

- marker: marker name
- genetic: genetic position estimated with R/qtl est.map 
- chrom: chromosome
- pos: physical position
- cross: the RIAILs supanels used to constuct the map (see Table S2)
- mapfunction = the map function used in the R/qtl est.map function
- intercross: intercross type (or pooled = consensus) (see Table S2)
- rec: rec-1 allele
- domain_type: domain inferred from Rockman & Kruglyak, 2009 https://doi.org/10.1371/journal.pgen.1000419
- marker_background: shared if segregating in both intercross, unshared if not

###
reference_comparison.R: Code to plot figure S2 + to calculate correlation with reference maps