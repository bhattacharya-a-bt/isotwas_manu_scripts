### This script provides sample code to simulate gene and isoform
### expression and quantitative phenotypes with different genetic architectures
### and uses isoTWAS and TWAS to map trait associations. The phenotype is simulated
### under Scenario 1 or 2 from Figure 3c of Bhattacharya et al 2022: (1)
### the simulated gene has an effect on the trait, but none of the observed
### isoforms have an effect, and (2) only one isoform of a gene has an effect
### on the trait with varying usage of the effect isoform (.1,.3,.5,.7,.9).

power.sim <- function(snps,
                      p.causal,
                      eqtl_h2,
                      het,
                      cor.struc,
                      n.tx,
                      gwas_arch,
                      gwas_var){
  snpsMat = as.matrix(snps$genotype[])
  snpsMat = apply(snpsMat,1,function(x) (x-mean(x))/sd(x))

  ### Simulate cis expression
  Z_qtl = simGeno(t(snps$genotypes[]), 400)
  p = ncol(Z_qtl)
  A_qtl = (t(Z_qtl) %*% Z_qtl) / (p-1)
  n = nrow(Z_qtl)
  p = ncol(Z_qtl)

  # GRM AND LD
  A = (Z_qtl %*% t(Z_qtl)) / (p-1)

  if (het == 'same'){
    b = simBeta(p.causal,
                eqtl_h2,
                n_snps = ncol(Z_qtl))
    b_qtls = matrix(b,
                    length(b),
                    n.tx)
  }
  if (het == 'different'){
    b_qtls = replicate(n.tx,simBeta(p.causal,
                                    eqtl_h2,
                                    n_snps = ncol(Z_qtl)))
  }

  # sim gene expression
  gexpr = apply(b_qtls,
                2,
                function(x) simTrait(Z_qtl %*% x, eqtl_h2))

  cor.mat = matrix(nrow = n.tx,
                   ncol = n.tx)

  if (cor.struc == 'low'){
    cor.mat = cov2cor(clusterGeneration::genPositiveDefMat(n.tx,
                                                           covMethod='c-vine',
                                                           eta=1000)$Sigma)
  } else {
    cor.mat = cov2cor(clusterGeneration::genPositiveDefMat(n.tx,
                                                           covMethod='c-vine',
                                                           eta=1)$Sigma)

  }

  gexpr.error = MASS::mvrnorm(400,mu = rep(0,n.tx),
                              Sigma = cor.mat)
  gexpr.tot = gexpr + gexpr.error
  colnames(gexpr.tot) = paste0('Isoform',1:ncol(gexpr))
  rownames(gexpr.tot) = rownames(Z_qtl) = paste0('s',1:nrow(gexpr))
  colnames(Z_qtl) = paste0('SNP',1:ncol(Z_qtl))
  id =  paste0('s',1:nrow(gexpr.tot))



  ### Run isoTWAS
  require(isoTWAS)
  model.iso = compute_isotwas(X = Z_qtl,
                              Y = gexpr.tot,
                              Y.rep = gexpr.tot,
                              R = 1,
                              id = as.factor(id),
                              omega_est = 'mean',
                              omega_nlambda = 10,
                              method = c('multi_enet',
                                         'univariate',
                                         'mvsusie',
                                         'mrce_lasso'),
                              predict_nlambda = 20,
                              family = 'gaussian',
                              scale = F,
                              alpha = 0.5,
                              nfolds = 5,
                              verbose = T,
                              par = F,
                              n.cores = 1,
                              tx_names = NULL,
                              seed = 1218,
                              run_all = F,
                              return_all = T,
                              tol.in = 1e-3,
                              maxit.in = 1e3,
                              coverage = .9)
  gexpr = as.matrix(rowSums(gexpr.tot))
  colnames(gexpr) = 'gene'

  uni_blup = univariate_blup(X = Z_qtl,
                             Y = gexpr,
                             Omega = 1,
                             scale = F,
                             alpha = .5,
                             nfolds = 5,
                             verbose = T,
                             par = F,
                             n.cores = 1,
                             tx_names = NULL,
                             seed = 1218)
  uni_enet = univariate_elasticnet(X = Z_qtl,
                                   Y = gexpr,
                                   Omega = 1,
                                   scale = F,
                                   alpha = .5,
                                   nfolds = 5,
                                   verbose = T,
                                   par = F,
                                   n.cores = 1,
                                   tx_names = NULL,
                                   seed = 1218)
  uni_susie = univariate_susie(X = Z_qtl,
                               Y = gexpr,
                               Omega = 1,
                               scale = F,
                               alpha = .5,
                               nfolds = 5,
                               verbose = T,
                               par = F,
                               n.cores = 1,
                               tx_names = NULL,
                               seed = 1218)



  univariate = list(uni_blup,uni_enet,uni_susie)
  rrr = unlist(sapply(univariate,
                      function(x) x[[1]]$R2))
  model.gene = univariate[[which(rrr == sample(c(max(rrr),
                                                 median(rrr)),
                                               1))]]


  #### SIM GWAS
  if (gwas_arch == 'gene_not_iso'){


    Z_gwas = simGeno(t(snps$genotypes[]), 20000)
    colnames(Z_gwas) = colnames(Z_qtl)
    gexpr_gwas = apply(b_qtls,
                       2,
                       function(x) simTrait(Z_gwas %*% x, eqtl_h2))
    gexpr_gwas_tot = as.numeric(rowSums(gexpr_gwas))

    if (gwas_var > 0){
      alpha = rnorm(1)
    }  else {alpha = 0}

    y = simTrait(gexpr_gwas_tot * alpha, gwas_var)
    gwas = regress(Z_gwas,y)

    gene_level_grex = as.numeric(model.gene[[1]]$Model$Weight %*%
                                   t(Z_gwas[,model.gene[[1]]$Model$SNP]))

    iso.mat = matrix(ncol = n.tx,
                     nrow = ncol(Z_qtl))
    rownames(iso.mat) = colnames(Z_qtl)
    for (col in 1:ncol(iso.mat)){
      mmm = model.iso[[1]][[col]]$Model
      for (row in 1:nrow(iso.mat)){
        if (rownames(iso.mat)[row] %in% mmm$SNP) {
          iso.mat[row,col] = mmm$Weight[mmm$SNP == rownames(iso.mat)[row]]
        }
      }
    }
    iso.mat[is.na(iso.mat)] = 0

    iso_level_grex = Z_gwas %*% iso.mat
    gene_p = coef(summary(lm(as.numeric(y) ~ as.numeric(gene_level_grex))))[2,4]
    iso_p = apply(iso_level_grex,2,
                  function(x){
                    coef(summary(lm(as.numeric(y) ~ x)))[2,4]
                  })
    iso_p = c(iso_p,
              coef(summary(lm(as.numeric(y) ~
                                as.numeric(rowSums(iso_level_grex)))))[2,4])
    acat_p = ACAT::ACAT(iso_p)
    stage_p = min(iso_p)
    if (gene_p <= 2.5e-6){
      if (any(c(acat_p,stage_p) <= 2.5e-6)){
        g_weights = rep(0,nrow(gwas))
        iso.mat = as.matrix(iso.mat[,which(iso_p[1:n.tx] <= 2.5e-6)])
        for (i in 1:length(g_weights)){
          if (paste0('SNP',i) %in% model.gene[[1]]$Model$SNP){
            g_weights[i] = model.gene[[1]]$Model$Weight[model.gene[[1]]$Model$SNP == paste0('SNP',i)]
          }
        }
        ct = conditional_test(g_weights,
                              iso.mat,
                              ld = A_qtl,
                              gwas$Beta/gwas$SE,
                              gene_name = 'Gene',
                              tx_name = paste0('Isoform',which(iso_p[1:n.tx] <=
                                                                 2.5e-6)))
      }
    }

  }

  if (gwas_arch == 'only_one_iso'){
    Z_gwas = simGeno(t(snps$genotypes[]), 20000)
    colnames(Z_gwas) = colnames(Z_qtl)
    props = c(.1,.3,.5,.7,.9)
    for (p in props){
      rrr = runif(n = n.tx-1,min = 0, max = 2)
      rrr = rrr/sum(rrr)
      rrr = rrr/(1/(1-p))
      prop.vec = c(p,rrr)

      gexpr_gwas = apply(b_qtls,
                         2,
                         function(x) simTrait(Z_gwas %*% x, eqtl_h2))
      gexpr_gwas_tot = drop(gexpr_gwas %*% prop.vec)

      if (gwas_var > 0){
        alpha = rnorm(1)
      }  else {alpha = 0}

      y = simTrait(gexpr_gwas[,1] * alpha,
                   gwas_var)
      gwas = regress(Z_gwas,y)
      gene_level_grex = as.numeric(model.gene[[1]]$Model$Weight %*%
                                     t(Z_gwas[,model.gene[[1]]$Model$SNP]))

      iso.mat = matrix(ncol = n.tx,
                       nrow = ncol(Z_qtl))
      rownames(iso.mat) = colnames(Z_qtl)
      for (col in 1:ncol(iso.mat)){
        mmm = model.iso[[1]][[col]]$Model
        for (row in 1:nrow(iso.mat)){
          if (rownames(iso.mat)[row] %in% mmm$SNP) {
            iso.mat[row,col] = mmm$Weight[mmm$SNP == rownames(iso.mat)[row]]
          }
        }
      }
      iso.mat[is.na(iso.mat)] = 0

      iso_level_grex = Z_gwas %*% iso.mat
      gene_p = coef(summary(lm(as.numeric(y) ~ as.numeric(gene_level_grex))))[2,4]
      iso_p = apply(iso_level_grex,2,
                    function(x){
                      coef(summary(lm(as.numeric(y) ~ x)))[2,4]
                    })
      iso_p = c(iso_p,
                coef(summary(lm(as.numeric(y) ~
                                  as.numeric(rowSums(iso_level_grex)))))[2,4])
      acat_p = ACAT::ACAT(iso_p)
      stage_p = min(iso_p)
      if (gene_p <= 2.5e-6){
        if (any(c(acat_p,stage_p) <= 2.5e-6)){
          g_weights = rep(0,nrow(gwas))
          iso.mat = as.matrix(iso.mat[,which(iso_p[1:n.tx] <= 2.5e-6)])
          for (i in 1:length(g_weights)){
            if (paste0('SNP',i) %in% model.gene[[1]]$Model$SNP){
              g_weights[i] = model.gene[[1]]$Model$Weight[model.gene[[1]]$Model$SNP == paste0('SNP',i)]
            }
          }

          ct = conditional_test(g_weights,
                                iso.mat,
                                ld = A_qtl,
                                gwas$Beta/gwas$SE,
                                gene_name = 'Gene',
                                tx_name = paste0('Isoform',which(iso_p[1:n.tx] <=
                                                                   2.5e-6)))
        }
      }

    }

  }

}
