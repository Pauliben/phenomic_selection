
#Name: Paul Adunola and Felipe Ferrao
#Date: Oct. 23, 2024
#Title: PHENOMIC-ASSISTED SELECTION: ASSESSMENT OF THE POTENTIAL OF NEAR-INFRARED SPECTROSCOPY FOR BLUEBERRY BREEDING

library(BGLR)
library(AGHmatrix)
library(dplyr)

########################################################

load("bb_nir_dat.rda",verbose = T)

########################################################

#Phenomic-Genomic Prediction models 

########################################################

#NIR - BayesB

nir_gs = function(trait_name,dat,reps){
  
  set.seed(123); seed=sample(1:1e4,reps)
  rep_r2_bb = list();
  
  #Select NIR data
  nir = dat[ , grepl( "X" , names(dat) ) ]
  #Scale NIR and remove null columns
  nir.scale = scale(nir)
  nir.scale = nir.scale[ , colSums(is.na(nir.scale))==0]
  y_vec = dat[,trait_name]
  dat.pred = data.frame(trait = y_vec,nir.scale);dat.pred = na.omit(dat.pred)
  X = as.matrix(dat.pred[,-1])
  y = dat.pred$trait

  nIter=6000; burnIn=1000
  
  for (i in 1:reps) {
    set.seed(seed[i])
    cat(paste(i," "))
    folds=sample(1:10,size=length(y),replace=T)
    y[is.nan(y)]<-NA
    r2_gp = vector();r2_bb = vector();r2_bl = vector()
    for (j in 1:max(folds)) {
      test = which(folds==j)
      yNA = y; yNA[test] = NA
      #Fit
      fmBB=BGLR(y=yNA,ETA=list( list(X=X,model='BayesB')),nIter=nIter,burnIn=burnIn,verbose=F)
      yhat = fmBB$yHat[test]
      r2_bb[j] = cor(yhat, y[test],use="complete")
    }
    rep_r2_bb[[i]] = r2_bb
  }
  return(rep_r2_bb)
}

trait_name = "Bloom"
nir_mod = nir_gs(trait_name,pheno_nir,5)
mean(unlist(nir_mod))

########################################################

#NIR + Pedigree - BayesB

multi_anir_gs = function(A,trait_name,dat,reps){
  
  #Select NIR data
  nir = dat[ , grepl( "X" , names(dat) ) ]
  #Scale NIR and remove null columns
  nir.scale = scale(nir)
  nir.scale = nir.scale[ , colSums(is.na(nir.scale))==0]
  y_vec = dat[,trait_name]
  dat.pred = data.frame(trait = y_vec,nir.scale); 
  rownames(dat.pred) = dat$genotype
  X = as.matrix(dat.pred[,-1])
  y = dat.pred$trait

  nIter=6000; burnIn=1000
  
  set.seed(123); seed=sample(1:1e4,reps)
  rep_r2_bb = list();
  ETA2=list(A=list(K=A,model='RKHS'),G=list(X=X,model='BayesB'))
  for (i in 1:reps) {
    set.seed(seed[i])
    cat(paste(i," "))
    folds=sample(1:10,size=length(y),replace=T)
    
    r2_bb = vector();
    for (j in 1:max(folds)) {
      test = which(folds==j);yNA = y; yNA[test] = NA
      # #Fit
      fit = BGLR(y=yNA,ETA=ETA2,nIter=nIter,burnIn=burnIn,verbose=F)
      yhat = fit$yHat[test]
      r2_bb[j] = cor(yhat, y[test],use="complete")
    }
    rep_r2_bb[[i]] = r2_bb;
  }
  return(rep_r2_bb)
}

trait_name = "Bloom"
anir_mod = multi_gnir_gs(ped_mat,trait_name,pheno_nir,5)
mean(unlist(anir_mod))

########################################################

#Genomic - BayesB

#impute NAs in SNP data
source("snp_imputation.R")
snp_imp = snp_imputation(snps,maf = 0.01,ploidy=4)

geno_gs = function(X,trait_name,dat,reps){
  
  set.seed(123); seed=sample(1:1e4,reps)
  rep_r2_bb = list();
  
  #Preparing data
  y_vec = dat[,trait_name]
  rownames(X) = NULL
  dat.pred = data.frame(trait = y_vec,X);
  y = dat.pred$trait
  
  nIter=6000; burnIn=1000
  
  for (i in 1:reps) {
    set.seed(seed[i])
    cat(paste(i," "))
    folds=sample(1:10,size=length(y),replace=T)
    y[is.nan(y)]<-NA
    r2_bb = vector();
    for (j in 1:max(folds)) {
      test = which(folds==j)
      yNA = y; yNA[test] = NA
      #Fit
      fmBB=BGLR(y=yNA,ETA=list( list(X=X,model='BayesB')),
                nIter=nIter,burnIn=burnIn,verbose=F)
      yhat = fmBB$yHat[test]
      r2_bb[j] = cor(yhat, y[test],use="complete")
    }
    rep_r2_bb[[i]] = r2_bb
  }
  return(rep_r2_bb)
}

trait_name = "Bloom"
geno_mod = geno_gs(snp_imp,trait_name,pheno_nir,5)
mean(unlist(geno_mod))

##################################################################

#NIR + Genomic - BayesB

geno_mat = Gmatrix(snps,ploidy = 4, maf = 0.01)

trait_name = "Bloom"
gnir_mod = multi_gnir_gs(geno_mat,trait_name,pheno_nir,5)
mean(unlist(gnir_mod))


