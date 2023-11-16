# The rec-1 RIAIL were genotyped by combining a 2b-Rad seq and a Massaray (MALDI-TOF) genotyping (see 
 Genetics article for details: https://doi.org/10.1101/2023.07.18.549456)

The raw 2b-Rad sequence and WGS parental sequences (.fq.gz) can be found on NCBI (PRJNA1037511 and PRJNA1040005, respectively). 

This directory contains:

- massaray_rec1riails.csv: the raw diploid genotype data for the Massaray (MALDI-TOF) genotyping (reference = N2 ws245)

- 2b-Rad_rec1riails.csv: The samples' information for the 2b-Rad sequence.
"RIAIL": the name of the RIAIL = a letter corresponding to the cross and a number. If 'replicate' is in the name, it means that there are technical replicates on the same line.
"lig.barcode": the barcode introduced during the ligation (see methods)
"index.barcode": the barcode introduced during the PCR step (see methods)
"lig": ligation barcode id (1 to 8)
"index": pcr barcode id (1 to 24)
"lib": the library number (1 or 2)
"code": a unique identifier corresponding to the name of the .fq.gz files
"cross": idenfier of the RIAIL subpanel differing from the initial cross (A to H; genetic background, rec-1 allele and direction of the cross vary; see methods) 
"parent1"
"parent2"

- trimming.R: R script used to verify the presence of the Bcg1 site in 2b-Rad sequence and trim them to the 36 bp inserted sample DNA fragment. The 2b-Rad .fq.gz files in NCBI (PRJNA1040005) must be trimmed to keep only the C. elegans DNA fragments.

- barcode directory contains the .txt files used to sort reads by RIAIL using stacks and a .xlsx file containing the sequence of the primers used
