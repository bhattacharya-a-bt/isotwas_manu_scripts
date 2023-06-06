library("optparse")
option_list <- list(
  make_option(c("-i", "--index"), action="store_true", default=TRUE,
              help="index",type='numeric'),
  make_option(c("-n", "--number"), action="store_true", default=TRUE,
              help="number of genes",type='numeric')
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

##ceiling = 1640
##number = 300

index = opt$index
number = opt$number

library(isotwas)
library(data.table)
library(bigsnpr)

start = (index-1)*number + 1
stop = (index)*number

gene_exp = fread('BigBrain_renormalized/gene_exp_bigbrain_042923_rint_resid.bed.gz',
                 data.table=F)
gene_exp = gene_exp[c(start:stop),]

tx_exp = fread('BigBrain_renormalized/tx_exp_bigbrain_042923_rint_resid.bed.gz',
               data.table = F)
tx_exp = subset(tx_exp,gid %in% gene_exp$gid)


dir.create('BigBrain_renormalized/tempfiles/',showWarnings = F,recursive = T)

snps = snp_attach('BigBrain_renormalized/BB_hmp_forQTL.rds')
isoform_fractions = 
  fread('BigBrain_0427/gene_summary_Ntx_IF.tsv')
isoform_fractions$gene_id = sapply(strsplit(isoform_fractions$gene_id,'[.]'),
                                   function(x) x[1])


runGTEx_isotwas = function(gene){
  
  tissue = 'Adult PFC'
  doneFile = paste0('BigBrain_renormalized/Models/Done/',
                    paste0('adultBB',gene,'.RDS'))
  
  gene_cur = subset(gene_exp,gid == gene)
  tx_cur = subset(tx_exp,gid == gene)
  
  chr = gene_cur$`#chr`[1]
  start = max(1,gene_cur$start[1] - 1e6)
  end = gene_cur$end[1] + 1e6
  
  gene_exp_mat = as.matrix(as.numeric(gene_cur[1,7:ncol(gene_cur)]))
  colnames(gene_exp_mat) = gene
  rownames(gene_exp_mat) = colnames(gene_cur)[7:ncol(gene_cur)]
  
  exp_mat = t(as.matrix(tx_cur[,7:ncol(gene_cur)]))
  colnames(exp_mat) = tx_cur$id
  rownames(exp_mat) = colnames(tx_cur)[7:ncol(tx_cur)]
  
  
  print(paste('TRAINING',gene,'FOR',tissue))
  
  if (!file.exists(doneFile)){
    
    
    isoform_fractions_this = subset(isoform_fractions,
                                    gene_id == gene)
    max_if = max(isoform_fractions_this$MaxIF)
    n_tx = ncol(exp_mat)
    print(paste0(gene,' has ',n_tx,' isoforms'))
    
    out_folder = 'BigBrain_renormalized/Models'
    out_folder_TWAS = file.path(out_folder,
                                'TWAS')
    out_folder_isoTWAS = file.path(out_folder,
                                   'isoTWAS')
    out_folder_comparison = file.path(out_folder,
                                      'Comparison')
    out_folder_done = file.path(out_folder,
                                'Done')
    out_folder_pred = file.path(out_folder,'Predictions')
    dir.create(out_folder,recursive = T,showWarnings = F)
    dir.create(out_folder_TWAS,recursive = T,showWarnings = F)
    dir.create(out_folder_isoTWAS,recursive = T,showWarnings = F)
    dir.create(out_folder_comparison,recursive = T,showWarnings = F)
    dir.create(out_folder_done,recursive = T,showWarnings = F)
    dir.create(out_folder_pred,recursive = T,showWarnings = F)
    
    snp_list = snps$map$marker.ID[snps$map$chromosome == chr &
                                    snps$map$physical.pos < end &
                                    snps$map$physical.pos > start]
    
    if (length(snp_list) > 1){
      bbfile = file.path('BigBrain_renormalized/tempfiles/',
                         paste0(index,gene,'_temp'))
      
      
      if (file.exists(paste0(bbfile,'.rds')) |
          file.exists(paste0(bbfile,'.bk'))){
        
        file.remove(c(paste0(bbfile,'.rds'),
                      paste0(bbfile,'.bk')))
        
      }
      
      snp_current = snp_attach(subset(snps,
                                      ind.col = which(snps$map$marker.ID %in%
                                                        snp_list),
                                      #ind.row = which(snps$fam$sample.ID %in% ids),
                                      backingfile = bbfile))
      
      ### TRAIN TWAS MODEL
      snpMat = as.matrix(snp_current$genotypes[])
      colnames(snpMat) = snp_current$map$marker.ID
      rownames(snpMat) = as.character(snp_current$fam$sample.ID)
      snpMat = snpMat[rownames(exp_mat),]
      gene_exp_mat = as.matrix(gene_exp_mat[rownames(exp_mat),])
      exp_mat = as.matrix(exp_mat[rownames(exp_mat),])
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
      twas_gene_level_prediction = last.mod[[1]]$Pred
      
      if (ncol(exp_mat) > 1){
        
        m.isot = compute_isotwas(X = snpMat,
                                 Y = exp_mat,
                                 Y.rep = exp_mat,
                                 gene_exp = as.numeric(gene_exp_mat),
                                 R = 1,
                                 id = rownames(snpMat),
                                 omega_est = 'replicates',
                                 omega_nlambda = 5,
                                 method = c('multi_enet',
                                            'univariate',
                                            'mrce_lasso',
                                            'spls',
                                            'joinet'),
                                 predict_nlambda = 5,
                                 family = 'gaussian',
                                 scale = F,
                                 alpha = 0.5,	
                                 nfolds = 5,
                                 verbose = F,
                                 par = T,
                                 n.cores = 4,
                                 tx_names = colnames(exp_mat),
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
        uni_v_multi_comp$joinet = as.numeric(unlist(uni_v_multi_comp$joinet))
        uni_v_multi_comp$spls = as.numeric(unlist(uni_v_multi_comp$spls))
        uni_v_multi_comp$multivariate = as.numeric(apply(uni_v_multi_comp,1,
                                                         function(x) max(as.numeric(x[c(2:5)]))))
        
        uni_v_multi_comp$Gene = gene
        uni_v_multi_comp$N_tx = n_tx
        uni_v_multi_comp$MaxIF = max_if
        
        fwrite(uni_v_multi_comp,
               file.path(out_folder_comparison,
                         paste0(gene,'_isoformcomparison.tsv.gz')))
        
        glmnet_r2_tot = tx2gene$R2[1]
        
        
      } else {
        
        Model = gene.model
        Model$Feature = colnames(exp_mat)
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
            aaa$Model$Feature = colnames(exp_mat)[j]
            aaa$Model = aaa$Model[,c('Feature','SNP','Chromosome','Position','Build',
                                     'ALT','REF','Weight','R2')]
            model.df = rbind(model.df,aaa$Model)
          }
          
          
        }
        
        
      }
      model.df = model.df[!is.na(model.df$R2),]
      model.df$Gene = gene
      
      colnames(tx_predictions) = colnames(exp_mat)
      rownames(tx_predictions) = rownames(exp_mat)
      pred_out = list(twas_gene_predictions = twas_gene_level_prediction,
                      isotwas_tx_predictions = tx_predictions,
                      tx2geneWeights = tx2gene$Weight_tx2gene)
      
      
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
      saveRDS(pred_out,file.path(out_folder_pred,paste0(gene,'.RDS')))
      fwrite(data.frame(Gene = gene,
                        R2_TWAS = unlist(gene.model$R2)[1],
                        R2_isoTWAS = ifelse(class(tx2gene) == 'character',
                                            unlist(gene.model$R2)[1],
                                            tx2gene$R2[1]),
                        R2_glmnet = glmnet_r2_tot,
                        N_tx = n_tx,
                        MaxIsoformFraction = max_if),
             file.path(out_folder_comparison,
                       paste0(gene,'_comparison.tsv.gz')))
      
      if (file.exists(paste0(bbfile,'.rds')) |
          file.exists(paste0(bbfile,'.bk'))){
        
        file.remove(c(paste0(bbfile,'.rds'),
                      paste0(bbfile,'.bk')))
        
      }
      
      file.create(doneFile)
      
    }
  }
}

tot = proc.time()
for (g in unique(gene_exp$gid)){
  a = proc.time()
  tryCatch(
    runGTEx_isotwas(g),
    error = function(e){
      print('Error, on to the next one')
    }
  )
  print(proc.time() - a)
}
proc.time() - tot
