### This script provides sample code to simulate gene and isoform
### expression and quantitative phenotypes with different genetic architectures
### and uses isoTWAS and TWAS to map trait associations. The phenotype is simulated
### under Scenario 3 from Figure 3c of Bhattacharya et al 2022: two effect isoforms
### for the same gene with varying ratios of effect sizes (-1 to 1).


#### create sims
require(bigsnpr)
snps = snp_attach(snp_readBed2('/path/to/chromosome1/bedfile',
                               backingfile = tempfile()))


### SNPs around CACNA1E
snps = snp_attach(subset(snps,
                         ind.col = which(snps$map$physical.pos > 181452653 - 1e6 &
                                           snps$map$physical.pos < 181777219 + 1e6)))


require(MOSTWAS)
require(isotwas)
p.causal = c(.001,0.01,.1)
eqtl_h2.loc = c(.1,.25)
cor.s = c('low','high')
n.tx = 2
het = c('same','different')
gwas_arch = 'isoform_alt'
var_explained = c(.1,.5)
params = expand.grid(p.causal,eqtl_h2.loc,het,cor.s,n.tx,gwas_arch,
                     var_explained)
colnames(params) = c('p.causal','eqtl.h2','heterogeneity',
                     'cor.s','n.tx','gwas_arch','gwas_var')

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
  Z_qtl = simGeno(t(snps$genotypes[]), 200)
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

  gexpr.error = MASS::mvrnorm(200,mu = rep(0,n.tx),
                              Sigma = cor.mat)
  gexpr.tot = gexpr + gexpr.error
  colnames(gexpr.tot) = paste0('Isoform',1:ncol(gexpr))
  rownames(gexpr.tot) = rownames(Z_qtl) = paste0('s',1:nrow(gexpr))
  colnames(Z_qtl) = paste0('SNP',1:ncol(Z_qtl))
  id =  paste0('s',1:nrow(gexpr.tot))



  ### Run isoTWAS
  require(isotwas)
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
  gexpr = as.matrix(log(rowSums(exp(gexpr.tot))))
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
  model.gene = univariate[[which(rrr == median(rrr))]]


  #### SIM GWAS

  Z_gwas = simGeno(t(snps$genotypes[]), 5000)
  colnames(Z_gwas) = colnames(Z_qtl)
  A_gwas = (Z_gwas %*% t(Z_gwas)) / (p-1)
  gexpr_gwas = apply(b_qtls,
                     2,
                     function(x) simTrait(Z_gwas %*% x, eqtl_h2))
  gexpr_gwas_tot = as.numeric(rowSums(gexpr_gwas))


  if (gwas_arch == 'isoform_alt'){
    alpha1 = rnorm(1)
    props = c(1,.5,.2,-.2,-.5,-1)
    for (p in props){
      alpha2 = p*alpha1
      y = simTrait(gexpr_gwas %*% c(alpha1,
                                    alpha2), gwas_var)
      #gwas = regress(Z_gwas,y)

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
      gene_p = ifelse(nrow(coef(summary(lm(as.numeric(y) ~ as.numeric(gene_level_grex))))) > 1,
                      coef(summary(lm(as.numeric(y) ~ as.numeric(gene_level_grex))))[2,4],
                      1)
      iso_p = apply(iso_level_grex,2,
                    function(x){
                      coef(summary(lm(as.numeric(y) ~ x)))[2,4]
                    })
      iso_p = c(iso_p,
                coef(summary(lm(as.numeric(y) ~
                                  as.numeric(rowSums(iso_level_grex)))))[2,4])

      full = as.data.frame(cbind(y,iso_level_grex))
      colnames(full)[1] = 'Gene'
      lm_all = lm(Gene~.,data=full)
      lm_null = lm(Gene~1,data=full)
      lrt_p = lmtest::lrtest(lm_all,lm_null)[2,5]


      acat_p = ACAT::ACAT(iso_p)
      stage_p = min(iso_p)
      df = data.frame(p.causal = p.causal,
                      eqtl_h2 = eqtl_h2,
                      het = het,
                      cor.struc = cor.struc,
                      n.tx = n.tx,
                      gwas_arch = gwas_arch,
                      gwas_var = gwas_var,
                      prop = p,
                      gene_p = gene_p,
                      acat_p = acat_p,
                      stage_p = stage_p,
                      chi_p = lrt_p)

      data.table::fwrite(df,
                         '/u/scratch/a/abtbhatt/isotwas_power_sim/ht_sim_res_altGWAS.tsv',
                         append=T,
                         row.names=F,
                         quote=F,
                         sep='\t')
    }

  }




}


power_sim_error = function(){
  tryCatch({power.sim(snps,
                      p.causal = params$p.causal[i],
                      eqtl_h2 = params$eqtl.h2[i],
                      het = params$heterogeneity[i],
                      cor.struc = params$cor.s[i],
                      n.tx = params$n.tx[i],
                      gwas_arch = params$gwas_arch[i],
                      gwas_var = params$gwas_var[i])},
           error = function(e){
             data.frame(p.causal = params$p.causal[i],
                        eqtl_h2 = params$eqtl.h2[i],
                        het = params$heterogeneity[i],
                        cor.struc = params$cor.s[i],
                        n.tx = params$n.tx[i],
                        gwas_arch = params$gwas_arch[i],
                        gwas_var = params$gwas_var[i],
                        prop = c(1,.5,.2,-.2,-.5,-1),
                        gene_p = NA,
                        acat_p = NA,
                        stage_p = NA
             )
           })}



a = proc.time()
ttt = pbapply::pbreplicate(1000,power_sim_error())
proc.time() - a
