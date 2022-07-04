### This script provides sample code to simulate gene and isoform
### expression with different genetic architectures and use gene-level
### and multivariate isoform-level prediction to compare predictive performance


#### create sims
require(bigsnpr)
snps = snp_attach(snp_readBed2('/path/to/chromosome1/bedfile',
                               backingfile = tempfile()))


### SNPs around CACNA1E
snps = snp_attach(subset(snps,
                         ind.col = which(snps$map$physical.pos > 181452653 - 1e6 &
                                           snps$map$physical.pos < 181777219 + 1e6)))


require(MOSTWAS)
require(isoTWAS)
n.qtl = c(100,500)
p.causal = c(0.001,0.01,.05,.1)
eqtl_h2.loc = c(.1,.25)
cor.s = c('low','high')
n.tx = c(2,5,10)
het = c('same','different')
params = expand.grid(n.qtl,p.causal,eqtl_h2.loc,het,cor.s,n.tx)
colnames(params) = c('n.qtl','p.causal','eqtl.h2','heterogeneity',
                     'cor.s','n.tx')

prediction.sim <- function(snps,n.qtl,p.causal,
                           eqtl_h2,het,
                           cor.struc,n.tx){
  snpsMat = as.matrix(snps$genotype[])
  snpsMat = apply(snpsMat,1,function(x) (x-mean(x))/sd(x))

  ### Simulate cis expression
  Z_qtl = simGeno(t(snps$genotypes[]), n.qtl)
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

  diag(cor.mat) = 1
  for (i in 1:(n.tx-1)){
    for (j in (i+1):n.tx){
      correl = ifelse(cor.struc == 'low',
                      runif(1,0,.1),
                      runif(1,.5,.8))
      correl = ifelse(rbinom(1,1,.5) == 1,
                      -correl,correl)
      cor.mat[i,j] = correl
      cor.mat[j,i] = correl
    }
  }

  gexpr.error = MASS::mvrnorm(n.qtl,mu = rep(0,n.tx),
                              Sigma = cor.mat)
  gexpr.tot = gexpr + gexpr.error
  colnames(gexpr.tot) = paste0('Isoform',1:ncol(gexpr))
  rownames(gexpr.tot) = rownames(Z_qtl) = paste0('s',1:nrow(gexpr))
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
                          omega_nlambda = 10,
                          method = c('mrce_lasso',
                                     'curds_whey',
                                     'multi_enet',
                                     'univariate',
                                     'finemap',
                                     'mvsusie'),
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
  mmm = as.matrix(model$R2[,-1])

  df = data.frame(n.qtl = n.qtl,
                  p.causal = p.causal,
                  eqtl.h2 = eqtl_h2,
                  het = het,
                  cor.structure = cor.struc,
                  n.tx = n.tx,
                  mrce_lasso = mean(unlist(mmm[,1])),
                  curds_whey = mean(unlist(mmm[,2])),
                  multi_enet = mean(unlist(mmm[,3])),
                  univariate = mean(unlist(mmm[,4])),
                  mvsusie = mean(unlist(mmm[,6])),
                  finemap = mean(unlist(mmm[,5])))


  return(df)

}

a = proc.time()
ttt = pbapply::pbreplicate(2,prediction.sim(snps,
                                            n.qtl = params$n.qtl[i],
                                            p.causal = params$p.causal[i],
                                            eqtl_h2 = params$eqtl.h2[i],
                                            het = params$heterogeneity[i],
                                            cor.struc = params$cor.s[i],
                                            n.tx = params$n.tx[i]))
proc.time() - a

df = as.data.frame(t(ttt[1,,]))
df$het = params$heterogeneity[i]
df$cor.structure = params$cor.s[i]
