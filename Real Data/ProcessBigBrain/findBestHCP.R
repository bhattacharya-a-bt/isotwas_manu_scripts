require(data.table)
setwd('/u/scratch/a/abtbhatt/renormalize_0429')

df = data.frame(k = seq(100,300,25),
                num_eGene = 0,
                num_totalGene = 0,
                num_eIsoform = 0,
                num_totalIsoform = 0)

for (i in 1:nrow(df)){
  print(i)
  k = df$k[i]
  gene = data.table::fread(paste0('hcp',
                      k
                      ,'_cis_gene.txt'))
  colnames(gene) = paste0('V',1:ncol(gene))
  df$num_totalGene[i] = length(unique(gene$V1))
  gene$FDR = p.adjust(gene$V12,'fdr')
  gene = subset(gene,FDR < 0.05)
  df$num_eGene[i] = length(unique(gene$V1))
  rm(gene)
  
  isoform = data.table::fread(paste0('hcp',
                         k,
                         '_cis_isoform.txt'))
  colnames(isoform) = paste0('V',1:ncol(isoform))
  df$num_totalIsoform[i] = length(unique(isoform$V1))
  isoform$FDR = p.adjust(isoform$V12,'fdr')
  isoform = subset(isoform,FDR < 0.05)
  df$num_eIsoform[i] = length(unique(isoform$V1))
  rm(isoform)
  
  
}

### FOR GENES
print(df$k[which.max(df$num_eGene)])
print(df$k[which.max(df$num_eIsoform)])

gene_bed = fread('gene_exp_bigbrain_042923_nohcp.bed.gz',data.table = F)
covs = fread('cov_hcp175_gene.txt')
require(limma)
gene_bed[,7:ncol(gene_bed)] = removeBatchEffect(gene_bed[,7:ncol(gene_bed)],
                                                covariates = t(covs[,-1]))
fwrite(gene_bed,
       'gene_exp_BB_050123_residualized.bed.gz',
       sep='\t',col.names=T,row.names=F)
rm(gene_bed)
rm(covs)

tx_bed = fread('tx_exp_bigbrain_042923_nohcp.bed.gz',data.table = F)
covs = fread('cov_hcp100_isoform.txt')
require(limma)
tx_bed[,7:ncol(tx_bed)] = removeBatchEffect(tx_bed[,7:ncol(tx_bed)],
                                            covariates = t(covs[,-1]))
fwrite(tx_bed,
       'tx_exp_BB_050123_residualized_isoform.bed.gz',
       sep='\t',col.names=T,row.names=F)
rm(tx_bed)
rm(covs)
