# Phenomic Selection in Blueberry

This folder contains the code and dataset used to compare phenomic and genomic selection models for predicting blueberry fruit quality traits.

There are three main files:

bb_nir_dat.rda: Contains the NIR, genomic, pedigree, and phenotypic data.
blueberry_phenomics.R: Contains the R code used to fit Bayesian BayesB, GBLUP, and GBLUP + BayesB models using BGLR-R with k-fold cross-validation.
snp_imputation.R: Contains the code for imputing missing SNP information.
