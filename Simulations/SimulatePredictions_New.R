prediction.sim <- function(snps,n.qtl,
                           p.causal,
                           prop_shared,
                           eqtl_h2,
                           eqtl_nongenetic_corr,
                           n.tx){
  
  snpsMat = as.matrix(snps$genotype[])
  snpsMat = apply(snpsMat,1,function(x) (x-mean(x))/sd(x))
  
  ### Simulate cis expression
  Z_qtl = simGeno(t(snps$genotypes[]), n.qtl)
  n = nrow(Z_qtl)
  p = ncol(Z_qtl)
  
  simB <- function(p.causal, 
                   prop_shared, 
                   eqtl_h2, 
                   n,
                   n.tx){
    
    B_matrix = matrix(nrow = n,
                      ncol = n.tx)
    
    n_qtls = max(1,floor(p.causal * n))
    which_shared = sample(1:n,floor(prop_shared * n_qtls))
    
    for (c in 1:ncol(B_matrix)){
      
      # select which SNPs are causal
      if (length(which_shared) > 0){
        c_qtls = sample((1:n)[-which_shared],
                        (n_qtls - length(which_shared))) 
      } else {
        c_qtls = sample(1:n,n_qtls) 
      }
      b_qtls = rep(0,n)
      
      # sample effects from normal prior
      b_qtls[c(c_qtls,which_shared)] = rnorm(n_qtls,
                                             mean = 0,
                                             sd = sqrt(eqtl_h2/n_qtls))
      B_matrix[,c] = b_qtls
      
    }
    
    return(B_matrix)
    
  }
  
  B_matrix = simB(p.causal, 
                  prop_shared, 
                  eqtl_h2, 
                  p,
                  n.tx)
  
  # sim gene expression with genetics
  gexpr_genetics = Z_qtl %*% B_matrix
  
  ### non genetics
  simulateCorr = function(dimension,low,high){
    
    mat = matrix(runif(dimension^2,low,high),
                 nrow = dimension,
                 ncol = dimension)
    diag(mat) = 1
    tmt = tcrossprod(mat)
    return(round(tmt/max(tmt),3))
    
  }
  
  V = simulateCorr(dimension = n.tx,
                   low = -.4,
                   high = .4)
  W = simulateCorr(dimension = n,
                   low = -.02,
                   high = .02)
  U = MixMatrix::rmatrixnorm(n=1,
                             mean = matrix(0,
                                           ncol = n.tx,
                                           nrow = n.qtl),
                             U = eqtl_nongenetic_corr*W,
                             V = eqtl_nongenetic_corr*V)
  
  E = t(mvnfast::rmvn(n = n.tx,
                      mu = rep(0,n.qtl),
                      sigma = (1 - eqtl_h2 - eqtl_nongenetic_corr) * 
                        diag(1,nrow = n.qtl)))
  
  gexpr.tot = gexpr_genetics + U + E
  colnames(gexpr.tot) = paste0('Isoform',1:ncol(gexpr.tot))
  rownames(gexpr.tot) = rownames(Z_qtl) = paste0('s',1:nrow(gexpr.tot))
  colnames(Z_qtl) = paste0('SNP',1:ncol(Z_qtl))
  id =  paste0('s',1:nrow(gexpr.tot))
  
  gene_expression = rowSums(gexpr.tot)
  
  
  
  ### Run isoTWAS
  require(isotwas)
  model = compute_isotwas(X = Z_qtl,
                          Y = gexpr.tot,
                          Y.rep = gexpr.tot,
                          R = 1,
                          id = as.factor(id),
                          omega_est = 'mean',
                          omega_nlambda = 5,
                          method = c('multi_enet',
                                     'mvsusie',
                                     'mrce_lasso',
                                     'univariate'),
                          predict_nlambda = 5,
                          family = 'gaussian',
                          scale = F,
                          alpha = 0.5,
                          nfolds = 5,
                          verbose = F,
                          par = F,
                          n.cores = 1,
                          tx_names = colnames(gexpr.tot),
                          seed = NULL,
                          run_all = F,
                          return_all = T,
                          tol.in = 1e-3,
                          maxit.in = 1e3,
                          coverage = .9)
  mmm = as.matrix(model$R2[,-1])
  isotwas_fitted = rep(0,n)
  for (i in 1:length(model$Model)){
    isotwas_fitted = isotwas_fitted + model$Model[[i]]$Pred
  }
  
  model_twas = glmnet::cv.glmnet(x = Z_qtl,
                                 y = gene_expression,
                                 nfolds = 5,
                                 keep = T,
                                 trace.it = T)
  twas_fitted = model_twas$fit.preval[,which.min(model_twas$cvm)]
  twas_r2 = summary(lm(gene_expression~twas_fitted))$adj.r.squared
  isotwas_r2 = cor(gene_expression,isotwas_fitted)^2
  
  
  df = data.frame(n.qtl = n.qtl,
                  p.causal = p.causal,
                  prop.shared = prop_shared,
                  eqtl.h2 = eqtl_h2,
                  eqtl_nongenetic_corr = eqtl_nongenetic_corr,
                  n.tx = n.tx,
                  mrce_lasso = mean(unlist(mmm[,3])),
                  multi_enet = mean(unlist(mmm[,1])),
                  univariate = mean(unlist(mmm[,4])),
                  mvsusie = mean(unlist(mmm[,2])),
                  isotwas_r2 = isotwas_r2,
                  twas_r2 = twas_r2)
  
  return(df)
  
}

#### create sims
require(bigsnpr)
require(MOSTWAS)
require(isotwas)
n.qtl = c(200,500,1000)[3:1]
p.causal = c(0.001,0.01,.05,.1)
prop.shared = c(0,.5,1)
eqtl_h2.loc = c(0.05,.1,.25)
eqtl_nongenetic_corr = c(.1,.25)
n.tx = c(2,5,10)[3:1]
params = expand.grid(n.qtl,
                     p.causal,
                     prop.shared,
                     eqtl_h2.loc,
                     eqtl_nongenetic_corr,
                     n.tx,
                     1:22)
colnames(params) = c('n.qtl',
                     'p.causal',
                     'prop.shared',
                     'eqtl.h2',
                     'heterogeneity',
                     'n.tx','chrom')

snp_folder = 'snp_files'
snp_files = file.path(snp_folder,
                      list.files(snp_folder))

chr = params$chrom[index]

done = data.table::fread('sims_pred_done.tsv')
done = subset(done,
              Chromosome == chr &
                n.qtl == params$n.qtl[index] &
                p.causal == params$p.causal[index] &
                prop.shared == params$prop.shared[index] &
                eqtl.h2 == params$eqtl.h2[index] &
                eqtl_nongenetic_corr == params$heterogeneity[index] &
                n.tx == params$n.tx[index])
if (nrow(done) > 1){print(paste0('Index ',index,'is done!'))}
if (nrow(done) == 0){
  
  bed_files = snp_files[grepl(paste0('_chr',chr,'_'),
                              snp_files)]
  bed_files = bed_files[grepl('.bed',
                              bed_files)]
  snps = snp_attach(snp_readBed2(bed_files[1],
                                 backingfile = tempfile()))
  a = proc.time()
  ttt = pbapply::pbreplicate(10,prediction.sim(snps,
                                               n.qtl = params$n.qtl[index],
                                               p.causal = params$p.causal[index],
                                               prop_shared = params$prop.shared[index],
                                               eqtl_h2 = params$eqtl.h2[index],
                                               eqtl_nongenetic_corr = 
                                                 params$heterogeneity[index],
                                               n.tx = params$n.tx[index]))
  proc.time() - a
  dir.create('sims_isotwas_prediction')
  df = as.data.frame(t(ttt[1,,]))
  df$PercDiff = ((unlist(df$isotwas_r2) - unlist(df$twas_r2))/abs(unlist(df$twas_r2))) * 100
  df$Chromosome = chr
  df$Gene = strsplit(strsplit(bed_files,'/')[[1]][8],'_')[[1]][1]
  data.table::fwrite(df,
                     'sims_res_isoTWAS_prediction.tsv',
                     append=T,
                     row.names=F,
                     quote=F,
                     sep='\t')
  data.table::fwrite(df[1,c('Chromosome',
                            'n.qtl','p.causal','prop.shared','eqtl.h2',
                            'eqtl_nongenetic_corr','n.tx')],
                     'sims_pred_done.tsv',
                     append=T,
                     row.names=F,
                     quote=F,
                     sep='\t')
  
}
