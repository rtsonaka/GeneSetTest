# GeneSetTest
Two-stage approach to test for gene, gene-set/pathway effects on longitudinal data


* Rcode_Two-Stage.R contains code to apply the two-stage approach as descibed in Tsonaka etal. (2012). 
* SupportFunctions.R contains all the supporting functions.
* Data.txt contains data simulated from a linear mixed effects model with several SNPs associated with progression.
* snps.txt contains simulated SNP data on 3 genes.

This code has been used in the simulation study descibed in Tsonaka et.al. (2012). It implements the two-stage approach and returns p-values with and without correcting for the variability of the empirical Bayes estimates for the random effects. 

Reference: Tsonaka, R., van der Helm - Mil, A., Houwing-Duistermaat, J. (2012). A two-stage mixed-effects model approach for gene-set analyses in candidate gene studies. Statistics in Medicine, 13, 1190 – 1202.

