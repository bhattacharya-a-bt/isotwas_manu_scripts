library("optparse")
option_list <- list(
  make_option(c("-i", "--index"), action="store_true", default=TRUE,
              help="index",type='numeric'),
  make_option(c("-n", "--number"), action="store_true", default=TRUE,
              help="number of genes",type='numeric')
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

##ceiling = 72
##number = 200


index = opt$index
number = opt$number

library(isotwas)
library(data.table)
library(bigsnpr)
library(vroom)
library(glmnet)
library(rrBLUP)
library(susieR)

library(limma)


exp_file="/residualized/leafcutter/estimates"
gen_file="/genotype/file/as/plink"
cov_file = '/covariate/file'
id_map = '/intron/to/gene/map'

load(id_map)

pheno = subset(pheno,ensemblID != '.')
pheno$ensemblID = sapply(strsplit(pheno$ensemblID,
                                  '[.]'),
                         function(x) x[1])

exp = vroom(exp_file)
exp = merge(pheno[,c('ID','ensemblID')],exp,by='ID')
rm(pheno)
gene = unique(exp$ensemblID)
start = (index-1)*number + 1
stop = (index)*number
gene_ids = gene[c(start:min(stop,length(gene)))]
exp = subset(exp,ensemblID %in% gene_ids)

### COVARIATE RESIDUALIZE
cov = vroom::vroom(cov_file)
cov_df = data.frame(t(cov[,-1]))
colnames(cov_df) = cov$id
cov_df = as.data.frame(apply(cov_df,2,as.numeric))
cov_df$sex = as.character(cov[cov$id == 'sex',-1])
mm = model.matrix(~.-1,data = cov_df)
rm(cov,cov_df)

exp[,6:ncol(exp)] = removeBatchEffect(as.matrix(exp[,6:ncol(exp)]),
                                      covariates = mm)



snps = snp_attach(snp_readBed2(paste0(gen_file,'.bed'),
                               backingfile = tempfile()))



runPEC_isotwas = function(gene){
  
    exp_current = exp[exp$ensemblID %in% gene,]
    
    if (nrow(exp_current) > 0){
      
      exp_mat = as.matrix(t(as.matrix(exp_current[,6:ncol(exp_current)])))
      colnames(exp_mat) = exp_current$ID
      rownames(exp_mat) = colnames(exp_current[6:ncol(exp_current)])
      chr = unique(exp_current$`#Chr`)[1]
      start = max(1,min(exp_current$start) - 1e6)
      end = max(exp_current$start) + 1e6
      snp_list = snps$map$marker.ID[snps$map$chromosome == chr &
                                      snps$map$physical.pos < end &
                                      snps$map$physical.pos > start]
      
      
      if (length(snp_list) > 1){
        snp_current = snp_attach(subset(snps,
                                        ind.col = which(snps$map$marker.ID %in%
                                                          snp_list)))
        
        ### Train isoTWAS model for LC
        exp_mat = as.matrix(exp_mat[snps$fam$sample.ID,])
        
        snpMat = as.matrix(snp_current$genotypes[])
        colnames(snpMat) = snp_current$map$marker.ID
        rownames(snpMat) = as.character(snp_current$fam$sample.ID)
        snpMat = snpMat[rownames(exp_mat),]
        
        if (ncol(exp_mat) > 1){
          
          m.isot = compute_isotwas(X = snpMat,
                                   Y = exp_mat,
                                   Y.rep = exp_mat,
                                   R = 1,
                                   gene_exp = NULL,
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
              aaa$Model = aaa$Model[,c('Feature','SNP','Chromosome','Position','Build',
                                       'ALT','REF','Weight','R2')]
              model.df = rbind(model.df,aaa$Model)
              
              
            }
          }
          
        } else {
          
          
          enet = univariate_elasticnet(X = snpMat,
                                       Y = exp_mat,
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
                                 Y = exp_mat,
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
                                   Y = exp_mat,
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
          gene.model$Feature = clu
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
          model.df = gene.model
          
        }
        
        model.df = model.df[!is.na(model.df$R2),]
        ld.cor = Matrix::Matrix(snp_cor(snp_current$genotypes))
        colnames(ld.cor) = rownames(ld.cor) = snp_current$map$marker.ID
        
        
        out_folder = '/out/folder/'
        out_folder_isoTWAS = '/isoTWAS_LC_HapMap3'
        out_folder_ld = '//LDMatrix'
        dir.create(out_folder,recursive = T,showWarnings = F)
        dir.create(out_folder_isoTWAS,recursive = T,showWarnings = F)
        dir.create(out_folder_ld,recursive = T,showWarnings = F)
        
        fwrite(model.df,file.path(out_folder_isoTWAS,
                                  paste0(gene,'_isoTWAS.tsv.gz')),
               sep='\t',
               col.names=T,
               row.names=F,
               quote=F)
        saveRDS(ld.cor,
                file.path(out_folder_ld,
                          paste0(gene,'_LDMatrix.RDS')))
      }
    }
}

doneFile = '/dummy/file/to/keep/track'
for (gene in gene_ids){
  done = fread(doneFile)
  print(gene)
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
