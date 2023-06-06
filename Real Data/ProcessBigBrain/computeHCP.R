### HCP 80 with seqPCs
hidden_convariate_linear <- function(F,Y,k,lambda,lambda2,lambda3,iter) {
  ## Use Example
  # hcp = hidden_convariate_linear(standardize(datSeq), standardize(t(datExpr)),k=10,iter = 100)
  #
  # function [Z,B,U,o,error,error1,error2,dz,db,du] = hidden_covariate_linear(F,Y,k,lambda,lambda2,lambda3,iter);
  # input:
  #      F: a matrix nxd of known covariates, where n is the number of
  #      subjects and d is the number of known covariates. *must be standardize (columns have 0 mean and constant SS).
  #      Y: a matrix of nxg of expression data (must be standardized (columns
  #      scaled to have constant SS and mean 0). ** use standardize function to standardize F and Y.
  #      k: number of inferred hidden components (k is an integer)
  #      lambda, lambda2, lambda3 are model parameters
  #      (optional) iter: number of iterations (default = 100);
  #
  #      note: k>0, lambda>0, lambda2>0, lambda3>0 must be set by the user based on the data at hand. one can set these values
  #      using cross-validation, by evaluating the "performance" of the  resulting residual data on a desired task.
  #      typically, if lambda>5, then hidden factors match the known covariates closely.
  #
  # objective:
  #
  # this function solves the following problem:
  # argmin_{Z,B,U}   ||Y-Z*B||_2 + \lambda*||Z-F*U||_2 + \lambda2*||B||_2 + \lambda_3||U||_2
  #
  # output:
  #      Z: matrix of hidden components, dimensionality: nxk
  #      B: matrix of effects of hidden components, dimensionality: kxg
  #      o: value of objective function on consecutive iterations.
  #
  # to use the residual data: Residual = Y - Z*B
  library(MASS)
  library(pracma)
  
  tol = 1e-6;
  
  U = matrix(0, nrow=dim(F)[2],k)
  Z = matrix(0, nrow=dim(F)[1],k)
  B = matrix(runif(dim(Z)[2]*dim(Y)[2]), nrow=dim(Z)[2], ncol=dim(Y)[2])
  F = as.matrix(F)
  
  n1 = dim(F)[1]
  d1 = dim(F)[2]
  
  n2 = dim(Y)[1]
  d2 = dim(Y)[2]
  
  if(n1!=n2)    stop("number of rows in F and Y must agree")
  
  if (k<1 | lambda<1e-6 | lambda2<1e-6 | lambda3<1e-6 ) {
    stop("lambda, lambda2, lambda3 must be positive and/or k must be an integer");
  }
  
  o = vector(length=iter)
  
  for (ii in 1:iter) {
    print(ii)
    o[ii] = sum((Y - Z%*%B)^2) + sum((Z -  F%*%U)^2)*lambda + 
      (sum(B^2))*lambda2 + lambda3*(sum(U^2));
    Z = (Y %*% t(B) + lambda * F %*%U) %*% ginv(B %*% t(B) + lambda * diag(dim(B)[1]))
    B = mldivide(t(Z) %*% Z + lambda2 * diag(dim(Z)[2]), (t(Z) %*% Y))
    U = mldivide(t(F) %*% F * lambda + lambda3 * diag(dim(U)[1]), lambda * t(F) %*% Z)
    
    if(ii > 1 &&  (abs(o[ii]-o[ii-1])/o[ii]) < tol)  break
  }
  
  error =  sum((Y - Z%*%B)^2) / sum(Y^2)  + sum((Z - F%*%U)^2)/sum((F%*%U)^2)
  error1 = sum((Y - Z%*%B)^2) / sum(Y^2);
  error2 = sum((Z - F%*%U)^2) / sum((F%*%U)^2);
  
  dz = Z%*%(B%*%t(B) + lambda*diag(dim(B)[1]))-(Y%*%t(B) + lambda*F%*%U);
  db = (t(Z)%*%Z + lambda2*diag(dim(Z)[2]))%*%B - t(Z)%*%Y;
  du = (t(F)%*%F*lambda + lambda3*diag(dim(U)[1]))%*%U-lambda*t(F)%*%Z;
  
  
  dataout = list(Z = Z, B = B, U = U)
  return(dataout)
}


standardize<- function(X)
{
  X = as.matrix(X)
  # n = dim(X)[1]
  # p = dim(X)[2]
  
  X = scale(X, center = TRUE, scale = F)
  # X = scale(X,center=FALSE, scale=sqrt(apply(X^2,2,sum)))
  
  # m = apply(X,2,mean)
  # st = sqrt(apply(X^2,2,sum));
  # st_mat = matrix(st, nrow = length(st), ncol = dim(X)[2], byrow=FALSE)
  # X2 = X / st_mat
  return (X)
}

require(data.table)
gene_bed = fread('/u/scratch/a/abtbhatt/renormalize_0429/gene_exp_bigbrain_042923_nohcp.bed.gz')
gene_cbt = as.matrix(gene_bed[,7:ncol(gene_bed)])

for (k in seq(100,300,by=25)){
  setwd('/u/scratch/a/abtbhatt/bigbrain_rnaseqqc/stats')
  file.remove('/u/scratch/a/abtbhatt/bigbrain_rnaseqqc/all_seq_stats.tsv')
  system('cat flipped* >> /u/scratch/a/abtbhatt/bigbrain_rnaseqqc/all_seq_stats.tsv')
  seqStats = fread('/u/scratch/a/abtbhatt/bigbrain_rnaseqqc/all_seq_stats.tsv')
  seqStats = subset(seqStats,
                    individualID != 'individualID')
  seqStatMat = as.matrix(seqStats[,-1])
  class(seqStatMat) = 'numeric'
  rownames(seqStatMat) = seqStats$individualID
  seqStatMat = seqStatMat[,which(apply(seqStatMat,2,var) != 0)]
  seqStatMat = seqStatMat[colnames(gene_bed)[7:ncol(gene_bed)],]
  
  
  outFolder = '/u/scratch/a/abtbhatt/renormalize_0429'
  hcp <- hidden_convariate_linear(standardize(seqStatMat), 
                                  standardize(t(gene_cbt)), 
                                  lambda=5,lambda2=1, lambda3=1, 
                                  k=k, iter=100)
  

  ### Aggregate covariates
  rnaseq = fread('/u/scratch/a/abtbhatt/renormalize_0429/rnaseq_050123_nodups.tsv')
  setwd('/u/project/gandalm/shared/GenomicDatasets-processed/PEC_AMP_AD_Arjun/Metadata')
  
  rnaseq = fread('rnaseq_20230501.csv')
  genotype = fread('genotype.csv')
  individual = fread('individual_20230421.csv')
  individual = subset(individual,individualID %in% 
                        rownames(seqStatMat))
  individual = individual[!duplicated(individual$individualID),]
  individual$reportedGender = tolower(individual$reportedGender)
  individual$sex_out = ifelse(individual$reportedGender != '',
                              individual$reportedGender,
                              individual$sex)
  individual$sex_out = as.factor(individual$sex_out)
  levels(individual$sex_out) = c('','F','F','M','M')
  
  covariates = data.frame(individualID = individual$individualID,
                          age = individual$ageDeath,
                          sex = individual$sex_out)
  covariates = covariates[match(rownames(seqStatMat),
                                covariates$individualID),]
  covariates$age[is.na(covariates$age)] = mean(covariates$age,
                                               na.rm=T)
  covariates$sex = ifelse(covariates$sex == 'F',0,1)
  covariates$age2 = covariates$age^2
  
  covariates = cbind(covariates,hcp$Z)
  colnames(covariates) = c('individualID','age','age2','sex',paste0('hcp',1:k))
  covariates = merge(covariates,individual[,c('individualID',
                                              'genotypingID')],by='individualID')
  
  
  gen_pc = fread('/u/project/gandalm/shared/GenomicDatasets-processed/PEC_AMP_AD_Arjun/Genotypes/all_pca.eigenvec')
  colnames(gen_pc)[2] = 'genotypingID'
  colnames(gen_pc)[3:22] = paste0('gPC',1:20)
  
  covariates = merge(covariates,gen_pc[,c('genotypingID',
                                          paste0('gPC',1:10))],
                     by='genotypingID')
  cov_out = data.frame(Covariate = colnames(covariates)[3:ncol(covariates)])
  cov_out = cbind(cov_out,t(covariates[,3:ncol(covariates)]))
  colnames(cov_out) = c('Covariate',covariates$individualID)
  cov_out = subset(cov_out,!Covariate %in% paste0('gPC',11:20))
  colnames(cov_out)[1] = 'id'
  fwrite(cov_out,
         paste0('/u/scratch/a/abtbhatt/renormalize_0429/cov_hcp',
                k,
                '_gene.txt'),
         sep='\t',
         col.names=T,
         row.names=F)
}




require(data.table)
tx_bed = fread('/u/scratch/a/abtbhatt/renormalize_0429/tx_exp_bigbrain_042923_nohcp.bed.gz')
tx_cbt = as.matrix(tx_bed[,7:ncol(tx_bed)])

for (k in seq(125,300,by=25)){
  setwd('/u/scratch/a/abtbhatt/bigbrain_rnaseqqc/stats')
  file.remove('/u/scratch/a/abtbhatt/bigbrain_rnaseqqc/all_seq_stats.tsv')
  system('cat flipped* >> /u/scratch/a/abtbhatt/bigbrain_rnaseqqc/all_seq_stats.tsv')
  seqStats = fread('/u/scratch/a/abtbhatt/bigbrain_rnaseqqc/all_seq_stats.tsv')
  seqStats = subset(seqStats,
                    individualID != 'individualID')
  seqStatMat = as.matrix(seqStats[,-1])
  class(seqStatMat) = 'numeric'
  rownames(seqStatMat) = seqStats$individualID
  seqStatMat = seqStatMat[,which(apply(seqStatMat,2,var) != 0)]
  seqStatMat = seqStatMat[colnames(tx_bed)[7:ncol(tx_bed)],]
  
  
  outFolder = '/u/scratch/a/abtbhatt/renormalize_0429'
  hcp <- hidden_convariate_linear(standardize(seqStatMat), 
                                  standardize(t(tx_cbt)), 
                                  lambda=5,lambda2=1, lambda3=1, 
                                  k=k, iter=100)
  
  
  ### Aggregate covariates
  rnaseq = fread('/u/scratch/a/abtbhatt/renormalize_0429/rnaseq_050123_nodups.tsv')
  setwd('/u/project/gandalm/shared/GenomicDatasets-processed/PEC_AMP_AD_Arjun/Metadata')
  
  genotype = fread('genotype.csv')
  individual = fread('individual_20230421.csv')
  individual = subset(individual,individualID %in% 
                        rownames(seqStatMat))
  individual = individual[!duplicated(individual$individualID),]
  individual$reportedGender = tolower(individual$reportedGender)
  individual$sex_out = ifelse(individual$reportedGender != '',
                              individual$reportedGender,
                              individual$sex)
  individual$sex_out = as.factor(individual$sex_out)
  levels(individual$sex_out) = c('','F','F','M','M')
  
  covariates = data.frame(individualID = individual$individualID,
                          age = individual$ageDeath,
                          sex = individual$sex_out)
  covariates = covariates[match(rownames(seqStatMat),
                                covariates$individualID),]
  covariates$age[is.na(covariates$age)] = mean(covariates$age,
                                               na.rm=T)
  covariates$sex = ifelse(covariates$sex == 'F',0,1)
  covariates$age2 = covariates$age^2
  
  covariates = cbind(covariates,hcp$Z)
  colnames(covariates) = c('individualID','age','age2','sex',paste0('hcp',1:k))
  covariates = merge(covariates,individual[,c('individualID',
                                              'genotypingID')],by='individualID')
  
  
  gen_pc = fread('/u/project/gandalm/shared/GenomicDatasets-processed/PEC_AMP_AD_Arjun/Genotypes/all_pca.eigenvec')
  colnames(gen_pc)[2] = 'genotypingID'
  colnames(gen_pc)[3:22] = paste0('gPC',1:20)
  
  covariates = merge(covariates,gen_pc[,c('genotypingID',
                                          paste0('gPC',1:10))],
                     by='genotypingID')
  cov_out = data.frame(Covariate = colnames(covariates)[3:ncol(covariates)])
  cov_out = cbind(cov_out,t(covariates[,3:ncol(covariates)]))
  colnames(cov_out) = c('Covariate',covariates$individualID)
  cov_out = subset(cov_out,!Covariate %in% paste0('gPC',11:20))
  colnames(cov_out)[1] = 'id'
  fwrite(cov_out,
         paste0('/u/scratch/a/abtbhatt/renormalize_0429/cov_hcp',
                k,
                '_isoform.txt'),
         sep='\t',
         col.names=T,
         row.names=F)
}



