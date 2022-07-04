### This script takes in a file that contains a matrix of transcript-isoform
### expression and builds an isoTWAS model. This script is built to be used
### in a job array through SGE or SLURM.

library("optparse")
option_list <- list(
  make_option(c("-f", "--filename"), action="store_true", default=TRUE,
              help="ENSG gene id",type='character')
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

filename = opt$filename[1]
gene = strsplit(strsplit(filename,'.RDS')[[1]][1],'/')[[1]][[9]]
tempDir = '/path/to/temp/folder'
dir.create(tempDir)
GCTA = 'gcta64'
plink = 'plink'
require(data.table)
markoff.file = '/path/to/file/keep/track/of/finished/genes'
markoff = fread(markoff.file,nThread = 1,fill=T)
markoff = subset(markoff,Gene == gene)
if (nrow(markoff) == 0){

  print(gene)

  require(data.table)
  exp_obj = readRDS(filename)

  eurID = rownames(exp_obj$Tx_Matrix)

  chr = strsplit(as.character(exp_obj$Coordinates$seqnames[1]),
                 'r')[[1]][2]
  if (chr %in% 1:22){
    start = max(c(1,min(exp_obj$Coordinates$start - 1e6)))
    end = max(exp_obj$Coordinates$end) + 1e6
    tx = colnames(exp_obj$Tx_Matrix)

    ### 1. EXTRACT PLINK FILE
    fwrite(data.frame(V1 = rownames(exp_obj$Tx_Matrix),
                      V2 = rownames(exp_obj$Tx_Matrix)),
           file.path(tempDir,paste0(gene,'_ids.txt')),
           sep='\t',
           col.names=F,
           row.names=F,
           quote=F,
           nThread = 1,
           na = 'NA')

    plink_start = paste0('/path/to/plink/prefix')
    system(paste0(plink,
                  ' --bfile ',plink_start,
                  ' --chr ',chr,
                  ' --from-bp ',start,
                  ' --to-bp ',end,
                  ' --keep ',file.path(tempDir,paste0(gene,'_ids.txt')),
                  ' --threads 1 --make-bed --out ',file.path(tempDir,gene)))

    if (file.exists(paste0(file.path(tempDir,gene),'.bed'))){

      cur.plink = paste0(file.path(tempDir,gene))
      final.id = fread(paste0(cur.plink,'.fam'),
                       nThread = 1)$V2
      exp_obj$Gene_Expression = as.matrix(rowSums(exp_obj$Tx_Matrix))
      colnames(exp_obj$Gene_Expression) = gene

      if (all(final.id %in% rownames(exp_obj$Tx_Matrix))){

        gene.exp.mat = log(as.matrix(exp_obj$Gene_Expression[final.id,])+1)
        colnames(gene.exp.mat) = gene

        exp.cov = as.data.frame(
          data.table::fread('/path/to/covariates',
                            header=T)
        )
        exp.cov = exp.cov[,c('id',rownames(gene.exp.mat))]
        exp.cov[6,-1] = as.numeric(ifelse(exp.cov[6,-1] == 'M',1,0))
        exp.cov.mat = as.matrix(exp.cov[,-1])
        class(exp.cov.mat) = 'numeric'

        gene.exp.mat = gene.exp.mat -
          t(exp.cov.mat) %*% MASS::ginv(exp.cov.mat %*%
                                          t(exp.cov.mat)) %*%
          exp.cov.mat %*% gene.exp.mat
        gene.exp.mat = apply(gene.exp.mat,2,function(x) {((x-mean(x))/sd(x))})

        tx.exp.mat = log(as.matrix(exp_obj$Tx_Matrix[final.id,])+1)
        colnames(tx.exp.mat) = colnames(exp_obj$Tx_Matrix)
        iso.cov = as.data.frame(
          data.table::fread('/path/to/covariates',
                            header=T)
        )
        iso.cov = iso.cov[,c('id',rownames(tx.exp.mat))]
        iso.cov[6,-1] = as.numeric(ifelse(iso.cov[6,-1] == 'M',1,0))
        iso.cov.mat = as.matrix(iso.cov[,-1])
        class(iso.cov.mat) = 'numeric'

        tx.exp.mat = tx.exp.mat -
          t(iso.cov.mat) %*% MASS::ginv(iso.cov.mat %*%
                                          t(iso.cov.mat)) %*%
          iso.cov.mat %*% tx.exp.mat
        tx.exp.mat = apply(tx.exp.mat,2,function(x) {((x-mean(x))/sd(x))})

        ### 2. HERITABILITY ANALYSES

        ### FIRST, GENE EXPRESSION

        system(paste0(GCTA,' --bfile ',cur.plink,
                      ' --thread-num 1 --autosome --make-grm-bin --out ',cur.plink))

        phen = as.data.frame(fread(paste0(cur.plink,'.fam'),
                                   nThread = 1))
        phen$V6 = as.numeric(gene.exp.mat[,1])
        phen = phen[,c(1,2,6)]
        fwrite(phen,paste0(cur.plink,'.phen'),sep='\t',col.names=F,
               row.names=F,quote=F)

        system(paste0(GCTA,
                      ' --reml --grm ',
                      cur.plink,
                      ' --pheno ',paste0(cur.plink,'.phen'),
                      ' --thread-num 1 --reml-no-constrain --reml-lrt 1 --out ',
                      cur.plink))

        hsq.df = data.frame(Name = c(exp_obj$Gene,
                                     as.character(exp_obj$Coordinates$tx_name)),
                            Category = c('Gene',
                                         rep('Tx',
                                             length(as.character(exp_obj$Coordinates$tx_name)))),
                            h2 = NA,
                            h2.se = NA,
                            h2.p = NA)

        if (file.exists(paste0(cur.plink,'.hsq'))){
          hsq = fread(paste0(cur.plink,'.hsq'),fill=T,nThread = 1)
          hsq.df$h2[1] = hsq$Variance[4]
          hsq.df$h2.se[1] = hsq$SE[4]
          hsq.df$h2.p[1] = hsq$Variance[9]
        } else {
          hsq.df$h2[1] = NA
          hsq.df$h2.se[1] = NA
          hsq.df$h2.p[1] = NA
        }


        ### SECOND, TX EXPRESSION
        for (i in 1:ncol(tx.exp.mat)){
          phen = as.data.frame(fread(paste0(cur.plink,'.fam'),
                                     nThread = 1))
          phen$V6 = as.numeric(tx.exp.mat[,i])
          phen = phen[,c(1,2,6)]
          fwrite(phen,paste0(cur.plink,'.phen'),sep='\t',col.names=F,
                 row.names=F,quote=F)

          if (file.exists(paste0(cur.plink,'.hsq'))){
            file.remove(paste0(cur.plink,'.hsq'))
          }

          system(paste0(GCTA,
                        ' --reml --grm ',
                        cur.plink,
                        ' --pheno ',paste0(cur.plink,'.phen'),
                        ' --thread-num 1 --reml-no-constrain --reml-lrt 1 --out ',
                        cur.plink))

          if (file.exists(paste0(cur.plink,'.hsq'))){
            hsq = fread(paste0(cur.plink,'.hsq'),nThread = 1,fill=T)
            hsq.df$h2[i+1] = hsq$Variance[4]
            hsq.df$h2.se[i+1] = hsq$SE[4]
            hsq.df$h2.p[i+1] = hsq$Variance[9]
          } else {
            hsq.df$h2[i+1] = NA
            hsq.df$h2.se[i+1] = NA
            hsq.df$h2.p[i+1] = NA
          }
        }

        hsq.df$h2 = abs(hsq.df$h2)
        ### 3. FIT MODEL FOR GENE EXPRESSION
        comp.gene = ifelse(is.na(hsq.df$h2.p[1]),1,hsq.df$h2.p[1])
        if (comp.gene <= 1){

          require(bigsnpr)
          require(isotwas)
          require(glmnet)
          require(rrBLUP)
          require(susieR)
          snps = snp_attach(snp_readBed2(paste0(cur.plink,'.bed'),
                                         backingfile = tempfile()))
          snpMat = snps$genotypes[]
          rownames(snpMat) = snps$fam$family.ID
          colnames(snpMat) = snps$map$marker.ID

          for (i in 1:ncol(snpMat)){

            mmm = snpMat[,i]
            mmm[which(is.na(mmm))] = mean(mmm,na.rm=T)
            snpMat[,i] = mmm

          }

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
                                       seed = 1218)

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
                                 seed = 1218)

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
                                   seed = 1218)

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

          colnames(last.mod[[1]]$Model)[1] = 'marker.ID'
          last.mod[[1]]$Model = merge(last.mod[[1]]$Model,
                                      snps$map[,c('chromosome',
                                                  'marker.ID',
                                                  'physical.pos',
                                                  'allele1',
                                                  'allele2')],
                                      by = 'marker.ID')
          colnames(last.mod[[1]]$Model) = c('SNP',
                                            'Weight',
                                            'Chromosome',
                                            'Position',
                                            'Allele1',
                                            'Allele2')
          last.mod[[1]]$Model =
            last.mod[[1]]$Model[order(last.mod[[1]]$Model$Position),]
          gene.model = list(Gene = gene,
                            Model = as.data.frame(last.mod[[1]]$Model),
                            R2 = max(r2),
                            R2.P = as.numeric(last.mod[[1]]$P))

          if (nrow(gene.model$Model) == 0){
            last.mod = susie
            colnames(last.mod[[1]]$Model)[1] = 'marker.ID'
            last.mod[[1]]$Model = merge(last.mod[[1]]$Model,
                                        snps$map[,c('chromosome',
                                                    'marker.ID',
                                                    'physical.pos',
                                                    'allele1',
                                                    'allele2')],
                                        by = 'marker.ID')
            colnames(last.mod[[1]]$Model) = c('SNP',
                                              'Weight',
                                              'Chromosome',
                                              'Position',
                                              'Allele1',
                                              'Allele2')
            last.mod[[1]]$Model =
              last.mod[[1]]$Model[order(last.mod[[1]]$Model$Position),]
            gene.model$Model = as.data.frame(last.mod[[1]]$Model)
          }

        } else {
          gene.model = NA
        }

        ## 4. ISOTWAS MODEL
        tx.h2 = subset(hsq.df,Category == 'Tx' & !is.na(h2.p) & h2 > 0)
        if (nrow(tx.h2) > 0 & any(tx.h2$h2.p < 0.05)){
          tx.h2 = subset(tx.h2,h2.p < .5)
          tx.exp.mat = as.matrix(tx.exp.mat[,tx.h2$Name])
          if (ncol(tx.exp.mat) > 1){

            require(bigsnpr)
            require(isotwas)
            require(glmnet)
            require(rrBLUP)
            require(susieR)
            snps = snp_attach(snp_readBed2(paste0(cur.plink,'.bed'),
                                           backingfile = tempfile()))
            snpMat = snps$genotypes[]
            rownames(snpMat) = snps$fam$sample.ID
            colnames(snpMat) = snps$map$marker.ID

            for (i in 1:ncol(snpMat)){

              mmm = snpMat[,i]
              mmm[which(is.na(mmm))] = mean(mmm,na.rm=T)
              snpMat[,i] = mmm

            }

            m.isot = compute_isotwas(X = snpMat,
                                     Y = tx.exp.mat,
                                     Y.rep = tx.exp.mat,
                                     R = 1,
                                     id = rownames(snpMat),
                                     omega_est = 'replicates',
                                     omega_nlambda = 10,
                                     method = c('mrce_lasso',
                                                'multi_enet',
                                                'univariate',
                                                'mvsusie'),
                                     predict_nlambda = 20,
                                     family = 'gaussian',
                                     scale = F,
                                     alpha = 0.5,
                                     nfolds = 5,
                                     verbose = F,
                                     par = F,
                                     n.cores = NULL,
                                     tx_names = colnames(tx.exp.mat),
                                     seed = 1218,
                                     run_all = F,
                                     return_all = T,
                                     tol.in = 1e-3,
                                     coverage = .9)
          } else {

            require(bigsnpr)
            require(isotwas)
            require(glmnet)
            require(rrBLUP)
            require(susieR)
            snps = snp_attach(snp_readBed2(paste0(cur.plink,'.bed'),
                                           backingfile = tempfile()))
            snpMat = snps$genotypes[]
            rownames(snpMat) = snps$fam$family.ID
            colnames(snpMat) = snps$map$marker.ID

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
                                         seed = 1218)

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
                                   seed = 1218)

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
                                     seed = 1218)

            r2 = unlist(c(enet[[1]]$R2,
                          blup[[1]]$R2,
                          susie[[1]]$R2))
            r2.min = which(r2 == max(r2))
            if (r2.min == 1){
              Model = enet
            }

            if (r2.min == 2){
              Model = blup
            }

            if (r2.min == 3){
              Model = susie
            }
            m.isot = list(Model = Model,
                          R2 = max(r2))

          }
          for (i in 1:length(m.isot$Model)){

            aaa = m.isot$Model[[i]]
            colnames(aaa$Model) = c('marker.ID','Weight')
            aaa$Model = merge(aaa$Model,
                              snps$map[,c('chromosome',
                                          'marker.ID',
                                          'physical.pos',
                                          'allele1',
                                          'allele2')],
                              by = 'marker.ID')
            colnames(aaa$Model) = c('SNP',
                                    'Weight',
                                    'Chromosome',
                                    'Position',
                                    'Allele1',
                                    'Allele2')
            aaa$Model = aaa$Model[order(aaa$Model$Position),]
            m.isot$Model[[i]] = aaa
            m.isot$Model[[1]]$Transcript = tx.h2$Name[1]

          }
        } else {
          m.isot = NA
        }


        final.list = list(Gene = gene,
                          Chromosome = chr,
                          Start = exp_obj$Start,
                          End = exp_obj$End,
                          Heritability = hsq.df,
                          Gene.Model = gene.model,
                          Isoform.Model = m.isot)
        outDir = '/path/to/out/folder'
        dir.create(outDir)
        saveRDS(final.list,file = file.path(outDir,paste0(gene,'.isotwas.wgt.RDS')))

      }
    }
  }

  fwrite(data.frame(Gene = gene,
                    Tissue = 'pe_fetal',
                    Done = 1),
         markoff.file,
         append=T,
         sep='\t',
         row.names=F,na=NA)
}

### DELETE TEMP FILES
setwd('/path/to/temp')
ddd = list.dirs()
ddd = ddd[grepl(gene,ddd)]
unlink(ddd,recursive=T)
