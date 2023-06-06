### This script takes in an 'index', defines a chunk of genes
### that is 'number' long and trains TWAS/isoTWAS models

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

exp_file="/path/to/normalized/isoform/expression"
gen_file="/path/to/snp/bed/file/" ### only include the prefix to .bed/.bim/.fam
hash_file='/path/to/metadata'
gtf_file='/path/to/gene/isoform/map'

library(isotwas)
library(data.table)
library(bigsnpr)

### Read in SNP file and restrict to HapMap3 SNPs
snps = snp_attach(snp_readBed2(paste0(gen_file,'.bed'),
                               backingfile = tempfile()))
hmp = fread('/path/to/hm.all.snp',
            header=F)
snps = snp_attach(subset(snps,ind.col = which(snps$map$marker.ID %in% hmp$V1)))
rm(hmp)

### Read in isoform expression data and restrict to the chunk of genes
exp = as.data.frame(fread(exp_file))
exp$ID = sapply(strsplit(exp$ID,'_'),
                function(x) x[1])
gtf_df <- fread(gtf_file)
gtf_df$Tx = sapply(strsplit(gtf_df$Tx,'_'),
                   function(x) x[1])
gtf_df = subset(gtf_df,Tx %in% exp$ID)
gene_ids = unique(gtf_df$Gene)
start = (index-1)*number + 1
stop = (index)*number
gene_ids = gene_ids[c(start:min(stop,length(gene_ids)))]
exp = subset(exp,ID %in% gtf_df$Tx[gtf_df$Gene %in% gene_ids])
gtf_df = subset(gtf_df,Gene %in% gene_ids)

### Read in gene expression data
gene_exp_data = fread('/path/to/covariate/normalized/gene/expression')
gene_exp_data = subset(gene_exp_data,TargetID %in% gene_ids)

runPEC_isotwas = function(gene){
  
  ### Subset to isoform of the gene
  gtf_df_current = gtf_df[gtf_df$Gene == gene,]
  gene_name = unique(gtf_df_current$Gene)[1]
  n_tx = length(gtf_df_current$Tx)
  exp_current = exp[exp$ID %in% gtf_df_current$Tx,]
  
  if (nrow(exp_current) > 0){
    
    ### Define the isoform expression matrix
    exp_mat = as.matrix(t(as.matrix(exp_current[,5:ncol(exp_current)])))
    colnames(exp_mat) = exp_current$ID
    rownames(exp_mat) = colnames(exp_current[5:ncol(exp_current)])
    
    ### Define the gene locus
    chr = unique(exp_current$`#Chr`)[1]
    start = max(1,min(exp_current$start) - 1e6)
    end = max(exp_current$start) + 1e6
    
    ### Subset to SNPs within the locus
    snp_list = snps$map$marker.ID[snps$map$chromosome == chr &
                                    snps$map$physical.pos < end &
                                    snps$map$physical.pos > start]
    
    
    if (length(snp_list) > 1){
      snp_current = snp_attach(subset(snps,
                                      ind.col = which(snps$map$marker.ID %in%
                                                        snp_list)))
      
      ### Train TWAS Model
      id_int = intersect(colnames(exp_current),
                         colnames(gene_exp_data))
      gene.exp.mat = t(as.matrix(gene_exp_data[gene_exp_data$TargetID == gene,
                                               5:ncol(gene_exp_data)]))
      gene.exp.mat = as.matrix(gene.exp.mat[id_int,])
      exp_mat = as.matrix(exp_mat[id_int,])
      
      snpMat = as.matrix(snp_current$genotypes[])
      colnames(snpMat) = snp_current$map$marker.ID
      rownames(snpMat) = as.character(snp_current$fam$sample.ID)
      snpMat = snpMat[rownames(exp_mat),]
      require(bigsnpr)
      require(isotwas)
      require(glmnet)
      require(rrBLUP)
      require(susieR)
      
      enet = univariate_elasticnet(X = snpMat,
                                   Y = gene.exp.mat,
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
                             Y = gene.exp.mat,
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
                               Y = gene.exp.mat,
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
      
      
      if (ncol(exp_mat) > 1){
        
        ### Train isoTWAS model
        m.isot = compute_isotwas(X = snpMat,
                                 Y = exp_mat,
                                 Y.rep = exp_mat,
                                 R = 1,
                                 gene_exp = as.numeric(gene.exp.mat),
                                 id = rownames(snpMat),
                                 omega_est = 'replicates',
                                 omega_nlambda = 5,
                                 method = c('mrce_lasso',
                                            'multi_enet',
                                            'univariate',
                                            'mvsusie'),
                                 predict_nlambda = 5,
                                 family = 'gaussian',
                                 scale = F,
                                 alpha = 0.5,	
                                 nfolds = 5,
                                 verbose = F,
                                 par = F,
                                 n.cores = NULL,
                                 tx_names = colnames(exp_mat),
                                 seed = 1218,
                                 run_all = F,
                                 return_all = T,
                                 tol.in = 1e-3,
                                 coverage = .9)
        
        tx2gene = m.isot$tx2gene_coef
        m.isot = m.isot$isotwas_mod
        
        ### Write out some information about the 
        uni_v_multi_comp = as.data.frame(m.isot$R2)
        uni_v_multi_comp$Transcript = 
          as.character(unlist(uni_v_multi_comp$Transcript))
        uni_v_multi_comp$mrce_lasso = 
          as.numeric(unlist(uni_v_multi_comp$mrce_lasso))
        uni_v_multi_comp$multi_enet = 
          as.numeric(unlist(uni_v_multi_comp$multi_enet))
        uni_v_multi_comp$univariate = 
          as.numeric(unlist(uni_v_multi_comp$univariate))
        uni_v_multi_comp$mvsusie = 
          as.numeric(unlist(uni_v_multi_comp$mvsusie))
        uni_v_multi_comp$multivariate = 
          apply(uni_v_multi_comp,1,
                function(x) max(as.numeric(x[c(2,3,5)])))
        
        uni_v_multi_comp$N_tx = n_tx
        uni_v_multi_comp$max_if = max(isoform_frac_cur$IsoformFraction)
        fwrite(uni_v_multi_comp,
               '/path/to/uni_v_multi_comp.tsv',
               sep='\t',
               append=T,
               row.names=F,
               quote=F)
        
        
      } else {
        
        ### If there is only 1 isoform, the gene and isoform model is the same
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
            this_model$Feature = exp_current$ID[1]
            model.df = rbind(model.df,this_model)
          }
          
        } else {
          
          colnames(aaa$Model) = c('marker.ID','Weight')
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
          aaa$Model = aaa$Model[,c('Feature','SNP','Chromosome',
                                   'Position','Build',
                                   'ALT','REF','Weight','R2')]
          model.df = rbind(model.df,aaa$Model)
          
          
        }
      }

      model.df = model.df[!is.na(model.df$R2),]
      ld.cor = Matrix::Matrix(snp_cor(snp_current$genotypes))
      colnames(ld.cor) = rownames(ld.cor) = snp_current$map$marker.ID
      
      
      out_folder = '/path/to/out/folder'
      out_folder_TWAS = '/TWAS_HapMap3'
      out_folder_isoTWAS = '/isoTWAS_HapMap3'
      out_folder_ld = '/LDMatrix'
      dir.create(out_folder,recursive = T,showWarnings = F)
      dir.create(out_folder_TWAS,recursive = T,showWarnings = F)
      dir.create(out_folder_isoTWAS,recursive = T,showWarnings = F)
      dir.create(out_folder_comparison,recursive = T,showWarnings = F)
      dir.create(out_folder_ld,recursive = T,showWarnings = F)
      
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
      saveRDS(ld.cor,
              file.path(out_folder_ld,
                        paste0(gene,'_LDMatrix.RDS')))
    }
  }
}

doneFile = '/dummy/file/to/keep/track'
for (gene in gene_ids){
  done = fread(doneFile)
  if (!(gene %in% done$Gene)){
    rm(done)
    tryCatch({
      runPEC_isotwas(gene)
    },error = function(e){print('Gene did not fit')})

    done_df = data.frame(Gene = gene,
                         Done = 1)
    fwrite(done_df,doneFile,
            append=T,row.names=F,
            sep='\t',quote=F)
  }
}
