library("optparse")
option_list <- list(
  make_option(c("-i", "--index"), action="store_true", default=TRUE,
              help="index",type='numeric'),
  make_option(c("-n", "--number"), action="store_true", default=TRUE,
              help="number of genes",type='numeric')
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


index = opt$index
number = opt$number

library(isotwas)
library(data.table)
library(bigsnpr)

base_folder = '/u/scratch/a/abtbhatt/GTEx_tx'
hash = fread(file.path(base_folder,
                       'file_location.tsv'))


start = (index-1)*number + 1
stop = (index)*number
hash = hash[start:stop,]

dir.create('/u/scratch/a/abtbhatt/gtextempfiles/',showWarnings = F,recursive = T)
tempfile = file.path('/u/scratch/a/abtbhatt/gtextempfiles',
                     paste0('allgtex_temp'))


gen_file = '/u/project/pasaniuc/pasaniucdata/GTEXv8_geno/GTEx.WGS.838.passOnly.geno0.05.hwe0.00001.dbsnp.SNPsOnly.NoAmbig.LDREF.bed'
snps = snp_attach(paste0(tempfile,'.rds'))
snps$fam$family.ID = stringr::str_replace_all(
  snps$fam$family.ID,
  '[.]',
  '-'
)

snps$fam$sample.ID = stringr::str_replace_all(
  snps$fam$sample.ID,
  '[.]',
  '-'
)

runGTEx_isotwas = function(i){
  
  gene = hash$Gene[i]
  tissue = hash$Tissue[i]
  
  isoform_fractions = fread(file.path('/u/scratch/a/abtbhatt/GTEx_tx',
                                      tissue,
                                      paste0(tissue,
                                             '_isoformfractions.tsv')))
  isoform_fractions = subset(isoform_fractions,
                             gene_id == gene)
  max_if = max(isoform_fractions$IsoformFraction)
  isoform_fractions = isoform_fractions[order(isoform_fractions$IsoformFraction,
                                              decreasing = T)]
  isoform_fractions$CumulativeSum = cumsum(isoform_fractions$IsoformFraction)
  
  gene_exp_obj = readRDS(file.path(base_folder,
                                   hash$File[i]))
  chr = hash$chromosome_name[i]
  gene_exp_mat = as.matrix(gene_exp_obj$Gene_Expression)
  tx_exp_mat = gene_exp_obj$Tx_Matrix
  n_tx = ncol(tx_exp_mat)
  
  if (ncol(tx_exp_mat) > 30){
    
    ccc = which( isoform_fractions$CumulativeSum <= .95 )
    if (length(ccc) == 0){
      ccc = c(1)
    } else {
      ccc = c(ccc,max(ccc) + 1)
    }
    
    keep_tx = isoform_fractions$isoform_id[ccc]
    tx_exp_mat = tx_exp_mat[,keep_tx]
    
  }
  
  
  if (nrow(gene_exp_obj$Coordinates) == ncol(tx_exp_mat) &
      is.null(colnames(tx_exp_mat))){
    colnames(tx_exp_mat) = unlist(gene_exp_obj$Coordinates$tx_name)
  }
  
  start = max(c(hash$start_position[i] - 1e6,1))
  end = hash$end_position[i] + 1e6
  
  snp_list = snps$map$marker.ID[snps$map$chromosome == chr &
                                  snps$map$physical.pos < end &
                                  snps$map$physical.pos > start]
  
  if (length(snp_list) > 1){
    bbfile = file.path('/u/scratch/a/abtbhatt/gtextempfiles/',
                       paste0(index,gene,tissue,'_temp'))
    
    
    if (file.exists(paste0(bbfile,'.rds')) |
        file.exists(paste0(bbfile,'.bk'))){
      
      file.remove(c(paste0(bbfile,'.rds'),
                    paste0(bbfile,'.bk')))
      
    }
    
    snp_current = snp_attach(subset(snps,
                                    ind.col = which(snps$map$marker.ID %in%
                                                      snp_list),
                                    backingfile = bbfile))
    
    ### TRAIN TWAS MODEL
    snpMat = as.matrix(snp_current$genotypes[])
    colnames(snpMat) = snp_current$map$marker.ID
    rownames(snpMat) = as.character(snp_current$fam$sample.ID)
    snpMat = snpMat[rownames(tx_exp_mat),]
    require(bigsnpr)
    require(isotwas)
    require(glmnet)
    require(rrBLUP)
    require(susieR)
    
    enet = univariate_elasticnet(X = snpMat,
                                 Y = gene_exp_mat,
                                 Omega = NA,
                                 family = 'gaussian',
                                 scale = F,
                                 alpha = 0.5,
                                 nfolds = 5,
                                 verbose = F,
                                 par = F,
                                 n.cores = NULL,
                                 tx_names = NULL,
                                 seed = NULL)
    
    blup = univariate_blup(X = snpMat,
                           Y = gene_exp_mat,
                           Omega = NA,
                           scale = F,
                           alpha = 0.5,
                           nfolds = 5,
                           verbose = F,
                           par = F,
                           n.cores = NULL,
                           tx_names = NULL,
                           seed = NULL)
    
    susie = univariate_susie(X = snpMat,
                             Y = gene_exp_mat,
                             Omega = NA,
                             scale = F,
                             alpha = 0.5,
                             nfolds = 5,
                             verbose = F,
                             par = F,
                             n.cores = NULL,
                             tx_names = NULL,
                             seed = NULL)
    
    r2 = unlist(c(enet[[1]]$R2,
                  blup[[1]]$R2,
                  susie[[1]]$R2))
    r2.min = which(r2 == max(r2))
    if (r2.min == 1){
      last.mod = enet
      if (nrow(enet[[1]]$Model) == 0){
        last.mod = susie
        r2.min = 3
      }
    }
    
    if (r2.min == 2){
      last.mod = blup
    }
    
    if (r2.min == 3){
      last.mod = susie
    }
    
    gene.model = last.mod[[1]]$Model
    gene.model$Feature = gene
    gene.model$R2 = last.mod[[1]]$R2
    gene.model$marker.ID = gene.model$SNP
    gene.model = merge(gene.model,snp_current$map,by='marker.ID')
    colnames(gene.model) = c('marker.ID','SNP',
                             'Weight','Feature','R2','Chromosome',
                             'genetic.dist','Position',
                             'ALT','REF')
    gene.model$Build = 'hg38'
    gene.model = gene.model[,c('Feature','SNP','Chromosome','Position','Build',
                               'ALT','REF','Weight','R2')]
    
    if (ncol(tx_exp_mat) > 1){
      
      
      m.isot = compute_isotwas(X = snpMat,
                               Y = tx_exp_mat,
                               Y.rep = tx_exp_mat,
                               gene_exp = as.numeric(gene_exp_mat),
                               R = 1,
                               id = rownames(snpMat),
                               omega_est = 'replicates',
                               omega_nlambda = 5,
                               method = c('multi_enet',
                                          'mvsusie',
                                          'univariate',
                                          'mrce_lasso'),
                               predict_nlambda = 5,
                               family = 'gaussian',
                               scale = F,
                               alpha = 0.5,	
                               nfolds = 5,
                               verbose = F,
                               par = T,
                               n.cores = 4,
                               tx_names = colnames(tx_exp_mat),
                               seed = NULL,
                               run_all = F,
                               return_all = T,
                               tol.in = 1e-3,
                               coverage = .9)
      tx2gene = m.isot$tx2gene_coef
      m.isot = m.isot$isotwas_mod
      uni_v_multi_comp = as.data.frame(m.isot$R2)
      uni_v_multi_comp$Transcript = as.character(unlist(uni_v_multi_comp$Transcript))
      uni_v_multi_comp$Tissue = tissue
      uni_v_multi_comp$mrce_lasso = as.numeric(unlist(uni_v_multi_comp$mrce_lasso))
      uni_v_multi_comp$multi_enet = as.numeric(unlist(uni_v_multi_comp$multi_enet))
      uni_v_multi_comp$univariate = as.numeric(unlist(uni_v_multi_comp$univariate))
      uni_v_multi_comp$mvsusie = as.numeric(unlist(uni_v_multi_comp$mvsusie))
      uni_v_multi_comp$multivariate = apply(uni_v_multi_comp,1,
                                            function(x) max(as.numeric(x[c(2,3,5)])))
      
      
      uni_v_multi_comp$Gene = gene
      uni_v_multi_comp$N_tx = ncol(as.matrix(gene_exp_obj$Tx_Matrix))
      fwrite(uni_v_multi_comp[,c('Gene','Transcript','N_tx',
                                 'mrce_lasso','multi_enet',
                                 'mvsusie','univariate',
                                 'Tissue','multivariate')],
             file.path(base_folder,
                       'uni_v_multi_comp_GTEx.tsv'),
             sep='\t',
             append=T,
             row.names=F,
             quote=F)
      
      
    } else {
      
      Model = gene.model
      Model$Feature = colnames(tx_exp_mat)
      m.isot = list(Model = list(Model = Model),
                    R2 = max(r2))
      tx2gene = 'Only 1 Tx'
      
    }
    
    model.df = data.frame(Feature = NA,
                          SNP = NA,
                          Chromosome = NA,
                          Position = NA,
                          Build = NA,
                          ALT = NA,
                          REF = NA,
                          Weight = NA,
                          R2 = NA)
    for (j in 1:length(m.isot$Model)){
      
      aaa = m.isot$Model[[j]]
      if (class(aaa) == 'data.frame'){
        this_model = aaa
        
        if (nrow(this_model) > 0){
          colnames(this_model)[colnames(this_model) == 'SNP'] = 
            c('marker.ID')
          this_model = merge(this_model,
                             snps$map,
                             by = 'marker.ID')
          this_model$R2 = unlist(this_model$R2)
          this_model = this_model[,c('Feature','marker.ID','Weight',
                                     'chromosome','physical.pos',
                                     'allele1','allele2','R2')]
          colnames(this_model) = c('Feature',
                                   'SNP',
                                   'Weight',
                                   'Chromosome',
                                   'Position',
                                   'REF',
                                   'ALT',
                                   'R2')
          this_model$Build = 'hg38'
          this_model = this_model[,c('Feature','SNP','Chromosome','Position','Build',
                                     'ALT','REF','Weight','R2')]
          model.df = rbind(model.df,this_model)
        }
        
      } else {
        
        colnames(aaa$Model) = c('marker.ID','Weight')
        if (nrow(aaa$Model) != 0){
          
          aaa$Model = merge(aaa$Model,
                            snps$map,
                            by = 'marker.ID')
          aaa$Model = aaa$Model[,c('marker.ID','Weight',
                                   'chromosome','physical.pos',
                                   'allele1','allele2')]
          colnames(aaa$Model) = c('SNP',
                                  'Weight',
                                  'Chromosome',
                                  'Position',
                                  'REF',
                                  'ALT')
          aaa$Model$R2 = aaa$R2
          aaa$Model$Build = 'hg38'
          aaa$Model$Feature = colnames(tx_exp_mat)[j]
          aaa$Model = aaa$Model[,c('Feature','SNP','Chromosome','Position','Build',
                                   'ALT','REF','Weight','R2')]
          model.df = rbind(model.df,aaa$Model)
        }
        
        
      }
      
      
    }
    model.df = model.df[!is.na(model.df$R2),]
    model.df$Gene = gene
    
    out_folder = file.path(base_folder,
                           'GTEx_Models_0423',
                           tissue)
    out_folder_TWAS = file.path(out_folder,
                                'TWAS_Models')
    out_folder_isoTWAS = file.path(out_folder,
                                   'isoTWAS_Models')
    out_folder_comparison = file.path(out_folder,
                                      'Comparison')
    dir.create(out_folder,recursive = T)
    dir.create(out_folder_TWAS,recursive = T)
    dir.create(out_folder_isoTWAS,recursive = T)
    dir.create(out_folder_comparison,recursive = T)
    
    fwrite(gene.model,file.path(out_folder_TWAS,
                                paste0(gene,'_TWAS.tsv.gz')),
           sep='\t',
           col.names=T,
           row.names=F,
           quote=F)
    saveRDS(list(SNP_to_Isoform = model.df,
                 Isoform_to_Gene = tx2gene),
            file.path(out_folder_isoTWAS,
                      paste0(gene,'_isoTWAS.RDS')))
    fwrite(data.frame(Gene = gene,
                      R2_TWAS = unlist(gene.model$R2)[1],
                      R2_isoTWAS = ifelse(class(tx2gene) == 'character',
                                          unlist(gene.model$R2)[1],
                                          tx2gene$R2[1]),
                      N_tx = n_tx,
                      MaxIsoformFraction = max_if),
           file.path(out_folder_comparison,
                     paste0(gene,'_comparison.tsv.gz')))
    
    if (file.exists(paste0(bbfile,'.rds')) |
        file.exists(paste0(bbfile,'.bk'))){
      
      file.remove(c(paste0(bbfile,'.rds'),
                    paste0(bbfile,'.bk')))
      
    }
    
  }
  
}

doneFile = file.path(base_folder,
                     'done_GTEx_models.tsv')

for (s in 1:nrow(hash)){
  
  done = vroom::vroom(doneFile)
  print(paste0('Gene: ',hash$Gene[s],
               ', Tissue: ',hash$Tissue[s]))
  if (!(hash$File[s] %in% done$File)){
    rm(done)
    tryCatch(
      runGTEx_isotwas(s),
      error = function(e){
        print('Error, on to the next one')
      }
    )
    
    done_df = data.frame(hash[s,c('Tissue','Gene','File')])
    fwrite(done_df,doneFile,append=T,row.names=F,sep='\t',quote=F)
  }
  
  
}
