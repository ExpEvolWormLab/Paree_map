# genotype

#RIAIL&P0_geno_ws245.csv

Filtered variants. The reference allele is the N2 allele (WS245). 

5 first columns are:
chromosome
position
reference allele
alternative allele
genotyping method (2b-Rad sequencing or MALDI-TOF massarray genotyping)

Then, the other columns are the genotype for each parental strain and RIAIL (0/0 | 1/1). 

The column name is the identifier of the RIAIL, and the first letter is the subpanel (see Table S2 for details)



#RIAL_geno_EEV1401polarized_qtlformat.csv

Genotype info in a format ready to be read with the R/qtl package in R.

The alleles are polarized toward the EEV1401 genetic background (BB).

The AA allele corresponds to the second parent's genetic background for the given cross (N2 or EEV1402)

Non-informative sites (monomorphic between parents of a given cross) are removed (NA value)

