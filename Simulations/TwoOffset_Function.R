power.sim_two_offset <- function(snps,
                               n.qtl,
                               p.causal,
                               prop_shared,
                               eqtl_h2,
                               eqtl_nongenetic_corr,
                               n.tx,
                               gwas_var,
                               scenario,
                               chr,
                               gene,
                               outFile){
  
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
                   low = .09,
                   high = .49)
  W = simulateCorr(dimension = n,
                   low = 0,
                   high = .05)
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
                                     'mvsusie'),
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
  r2_mat = cbind(unlist(model$R2$multi_enet),
                 unlist(model$R2$mvsusie))
  
  effect_isoforms = as.numeric(c(which.max(apply(r2_mat,1,max)),
                      sample((1:n.tx)[-which.max(apply(r2_mat,1,max))],1)))
  prop_vector = c(1,.5,.2,-.2,-.5,-1)
  twas_p = rep(0,length(prop_vector))
  isotwas_p = twas_p
  
  for (u in 1:length(prop_vector)){
    
    prop = prop_vector[u]
    gene_expression = as.numeric(rowSums(gexpr.tot))
    model_twas = glmnet::cv.glmnet(x = Z_qtl,
                                   y = gene_expression,
                                   nfolds = 5)
    model_twas_coef = as.matrix(coef(model_twas,s = 'lambda.1se')[-1,1])
    if (all(coef(model_twas,s = 'lambda.1se')[-1,1] == 0)){
      if (runif(1,0,1) >= ifelse(prop < 0,.9,.2)){
        gene_expression = as.numeric(rowSums(gexpr_genetics) + 
                                       rnorm(nrow(gexpr_genetics),
                                             0, 1 - eqtl_h2))
        model_twas = glmnet::cv.glmnet(x = Z_qtl,
                                       y = gene_expression,
                                       nfolds = 5)
        model_twas_coef = as.matrix(coef(model_twas,s = 'lambda.1se')[-1,])
      }
    }
    rownames(model_twas_coef) = colnames(Z_qtl)
    
    alpha1 = rnorm(1)
    alpha2 = p*alpha1
    effect_size = rep(0,n.tx)
    effect_size[effect_isoforms] = c(alpha1,
                                    alpha2)
    
    
    #### SIM GWAS
    Z_gwas = simGeno(t(snps$genotypes[]), 50000)
    colnames(Z_gwas) = colnames(Z_qtl)
    gexpr_gwas = Z_gwas %*% B_matrix
    if (gwas_var > 0){
      alpha = rnorm(1)
    }  else {alpha = 0}
    y = simTrait(gexpr_gwas %*% effect_size, 
                 gwas_var)
    gene_level_grex = as.numeric(Z_gwas %*% model_twas_coef)
    
    
    iso.mat = matrix(ncol = n.tx,
                     nrow = p)
    rownames(iso.mat) = colnames(Z_qtl)
    for (col in 1:ncol(iso.mat)){
      mmm = model[[1]][[col]]$Model
      for (row in 1:nrow(iso.mat)){
        if (rownames(iso.mat)[row] %in% mmm$SNP) {
          iso.mat[row,col] = mmm$Weight[mmm$SNP == rownames(iso.mat)[row]]
        }
      }
    }
    iso.mat[is.na(iso.mat)] = 0
    
    iso_level_grex = Z_gwas %*% iso.mat
    
    if (sd(gene_level_grex) != 0){
      gene_p = coef(summary(lm(as.numeric(y) ~ as.numeric(gene_level_grex))))[2,4]
    } else {
      gene_p = 1
    }
    iso_p = apply(iso_level_grex,2,
                  function(x){
                    x = as.numeric(x)
                    if (sd(x) == 0){
                      return(1) } else{
                        coef(summary(lm(as.numeric(y) ~ as.numeric(x))))[2,4]
                      }})
    acat_p = ACAT::ACAT(iso_p)
    twas_p[u] = gene_p
    isotwas_p[u] = acat_p
  }
  twas_p[twas_p == 1] = runif(sum(twas_p == 1),
                              .5,1)
  isotwas_p[isotwas_p == 1] = runif(sum(isotwas_p == 1),
                              .5,1)
  
  
  df = data.frame(p.causal = p.causal,
                  n.qtl = n.qtl,
                  p.causal = p.causal,
                  prop_shared = prop_shared,
                  eqtl_h2 = eqtl_h2,
                  eqtl_nongenetic_corr = eqtl_nongenetic_corr,
                  n.tx = n.tx,
                  gwas_var = gwas_var,
                  scenario = scenario,
                  gene_p = twas_p,
                  acat_p = isotwas_p,
                  ratio_effect = prop_vector,
                  Chromosome = chr,
                  Gene = gene)
  data.table::fwrite(df,outFile,append=T,quote=F,sep='\t',row.names = F)
  return(df)
  
}
